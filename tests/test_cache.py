"""验证完整入口优先使用已有数据，只有缺失文件才进入下载。"""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from pipeline.cache import available_stages, prepare_existing_data, seed_raw_cache
from pipeline.io import download, dump, sha
from pipeline.reuse import INPUTS

ROOT = Path(__file__).resolve().parents[1]


class CacheTests(unittest.TestCase):
    """用小文件检查完整计数复用、部分准备复用和原始文件核验。"""

    def setUp(self):
        """为每个用例生成隔离的新旧目录，不接触公开大文件。"""
        self.temp = tempfile.TemporaryDirectory()
        self.source = Path(self.temp.name)/'old'; self.source.mkdir()
        self.work = Path(self.temp.name)/'new'; self.work.mkdir()
        self.config = json.loads((ROOT/'config.json').read_text(encoding='utf-8'))
        dump(self.source/'run_config.json',self.config)
        dump(self.source/'owner.json',dict(source='official_public_downloads'))

    def tearDown(self):
        """清理当前测试生成的临时文件。"""
        self.temp.cleanup()

    def make_prepared(self, stages):
        """构造指定阶段的完整文件和摘要记录。"""
        for stage in stages:
            outputs = {}
            for name in INPUTS[stage]:
                p = self.source/name; p.parent.mkdir(parents=True,exist_ok=True)
                p.write_bytes(b'test data'); outputs[name] = sha(p)
            dump(self.source/'stages'/f'{stage}.json',dict(stage=stage,outputs=outputs))

    def test_prepared_never_downloads(self):
        """聚合计数齐全时不调用下载，也不扫描庞大的原始数据目录。"""
        self.make_prepared(INPUTS)
        with patch('urllib.request.urlopen',side_effect=AssertionError('不得联网')), \
             patch('pipeline.cache.seed_raw_cache',side_effect=AssertionError('不应重读原始文件')):
            found = prepare_existing_data(self.source,self.work,self.config)
        self.assertEqual(found,set(INPUTS))
        self.assertFalse((self.work/'results/mouse_results.csv.gz').exists())

    def test_partial_prepared(self):
        """只有人类输入准备完成时单独复用，不要求重新下载或聚合人类计数。"""
        self.make_prepared(['prepare_human'])
        found = prepare_existing_data(self.source,self.work,self.config)
        self.assertEqual(found,{'prepare_human'})
        self.assertTrue((self.work/'prepared/human/expression.csv.gz').is_file())
        self.assertFalse((self.work/'prepared/mouse_whole_counts.csv.gz').exists())

    def test_raw_skips_network(self):
        """完整原始文件导入后，下载函数直接返回，不访问远端。"""
        p = self.source/'raw/data.gz'; p.parent.mkdir(); p.write_bytes(b'complete')
        dump(p.with_name(p.name+'.download.json'),dict(url='https://example.invalid/data',
             sha256=sha(p),bytes=p.stat().st_size))
        seed_raw_cache(self.source,self.work)
        with patch('urllib.request.urlopen',side_effect=AssertionError('不应再次下载')):
            returned = download('https://example.invalid/data',self.work/'raw/data.gz')
        self.assertEqual(returned.read_bytes(),b'complete')

    def test_bad_raw_stops(self):
        """缓存损坏应停止而不是先覆盖旧文件重新下载。"""
        p = self.source/'raw/data.gz'; p.parent.mkdir(); p.write_bytes(b'changed')
        dump(p.with_name(p.name+'.download.json'),dict(url='x',sha256='wrong',bytes=7))
        with self.assertRaisesRegex(ValueError,'核验失败'):
            seed_raw_cache(self.source,self.work)

    def test_partial_is_not_linked(self):
        """断点文件复制后续写，不应改变原目录中的残片。"""
        p = self.source/'raw/data.part'; p.parent.mkdir(); p.write_bytes(b'abc')
        dump(p.with_name('data.partial.json'),dict(url='x',validator='tag'))
        seed_raw_cache(self.source,self.work)
        with (self.work/'raw/data.part').open('ab') as f:
            f.write(b'def')
        self.assertEqual(p.read_bytes(),b'abc')

    def test_unregistered_existing_raw(self):
        """有已知文件但没有来源记录时，不绕开它重复下载同一份数据。"""
        p = self.source/'raw/mouse.rds.gz'; p.parent.mkdir(); p.write_bytes(b'existing')
        with self.assertRaisesRegex(ValueError,'未重新下载'):
            seed_raw_cache(self.source,self.work)

    def test_resume_does_not_change_source(self):
        """恢复运行必须保持原数据来源。"""
        self.make_prepared(INPUTS)
        prepare_existing_data(self.source,self.work,self.config)
        another = Path(self.temp.name)/'other'; another.mkdir()
        with self.assertRaisesRegex(ValueError,'不能更换'):
            prepare_existing_data(another,self.work,self.config)


if __name__ == '__main__':
    unittest.main()
