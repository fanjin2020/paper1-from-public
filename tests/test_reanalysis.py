"""输入复用、模拟布局、按对比汇总和审阅包的回归测试。"""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile
import numpy as np
import pandas as pd
from pipeline.io import dump, sha
from pipeline.reuse import INPUTS, collect, import_inputs
from pipeline.simulation_design import synthetic, matrices, prepare, canonical_actual
from pipeline.audit import clinical_pairs
import run

ROOT = Path(__file__).resolve().parents[1]


class ReanalysisTests(unittest.TestCase):
    """不访问网络，用小型文件检查必须停止的分支和正常续跑。"""

    def setUp(self):
        """每个测试独立建立目录，避免共享缓存让失败用例意外通过。"""
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source = self.root / 'source'
        self.work = self.root / 'new'
        self.source.mkdir()
        self.work.mkdir()
        self.config = json.loads((ROOT / 'config.json').read_text(encoding='utf-8'))
        dump(self.source / 'owner.json', dict(source='official_public_downloads', signature='test'))
        dump(self.source / 'run_config.json', self.config)
        for stage, names in INPUTS.items():
            outputs = {}
            for name in names:
                p = self.source / name
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_bytes(b'small input\n')
                outputs[name] = sha(p)
            dump(self.source / 'stages' / f'{stage}.json', dict(stage=stage, outputs=outputs))

    def tearDown(self):
        """移除本测试创建的小型临时目录。"""
        self.temp.cleanup()

    def test_reuse_without_network(self):
        """输入导入不请求网络，来源文件哈希和原状态保持不变。"""
        before = {p.relative_to(self.source).as_posix(): sha(p) for p in self.source.rglob('*') if p.is_file()}
        with patch('urllib.request.urlopen', side_effect=AssertionError('不允许下载')):
            record = import_inputs(self.source, self.work, self.config)
        self.assertEqual(len(record['files']), 8)
        for name, digest in before.items():
            self.assertEqual(sha(self.source / name), digest)
        self.assertFalse((self.work / 'results' / 'mouse_results.csv.gz').exists())

    def test_reuse_resume(self):
        """相同输入的续跑不重新复制，也不需要来源模型结果。"""
        import_inputs(self.source, self.work, self.config)
        with patch('shutil.copy2', side_effect=AssertionError('不能再次复制')):
            import_inputs(self.source, self.work, self.config)

    def test_changed_input(self):
        """任何一个计数被改动，都必须在导入之前停止。"""
        (self.source / 'prepared/mouse_whole_counts.csv.gz').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'SHA-256'):
            collect(self.source, self.config)
        self.assertFalse((self.work / 'input_reuse.json').exists())

    def test_missing_receipt_entry(self):
        """文件存在但没有阶段身份记录，不能作为合格输入。"""
        dump(self.source / 'stages/prepare_mouse.json', dict(stage='prepare_mouse', outputs={}))
        with self.assertRaisesRegex(ValueError, '未登记'):
            collect(self.source, self.config)

    def test_parameter_change(self):
        """不同皮层定义不能复用同一聚合计数。"""
        altered = dict(self.config, cortex_labels=['another_region'])
        with self.assertRaisesRegex(ValueError, '参数不一致'):
            collect(self.source, altered)

    def test_model_change_allowed(self):
        """改变本次模型设置不会要求重下原始数据。"""
        record = collect(self.source, dict(self.config, ebayes_robust=True))
        self.assertEqual(len(record['files']), 8)

    def test_no_inplace_output(self):
        """来源目录与输出目录相同或嵌套时拒绝操作。"""
        for dst in [self.source, self.source / 'nested', self.root]:
            with self.assertRaises(ValueError):
                import_inputs(self.source, dst, self.config)

    def test_reuse_tampered_destination(self):
        """续跑时检查已导入文件，不能信任一次写过的状态。"""
        import_inputs(self.source, self.work, self.config)
        (self.work / 'results/panel.csv').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, '已改变'):
            import_inputs(self.source, self.work, self.config)

    def test_path_escape(self):
        """可选收据中的路径穿越也必须被拒绝。"""
        receipt = self.source / 'stages/prepare_human.json'
        record = json.loads(receipt.read_text(encoding='utf-8'))
        record['outputs']['prepared/human/../../../outside.csv'] = 'bad'
        dump(receipt, record)
        with self.assertRaisesRegex(ValueError, '越界'):
            collect(self.source, self.config)

    def test_clinical_vs_conditional(self):
        """临床状态和条件临床状态不能合并成一个显著集合。"""
        h = pd.DataFrame([dict(gene_symbol='A', context='Ast', spec=c, contrast=c, model=m,
                              beta=1., SE=.2, q_joint=.01) for c in
                          ['clinical_status', 'clinical_status_adjusted_for_braak'] for m in ['base', 'quality']])
        z = clinical_pairs(h)
        self.assertEqual(len(z), 1)
        self.assertEqual(z.contrast.iloc[0], 'clinical_status')

    def test_duplicate_endpoint(self):
        """重复检验位置不能被merge扩增后继续汇总。"""
        h = pd.DataFrame([dict(gene_symbol='A', context='Ast', spec='x', contrast='clinical_status', model=m)
                          for m in ['base', 'base', 'quality']])
        with self.assertRaises(pd.errors.MergeError):
            clinical_pairs(h)

    def test_exact_allocation(self):
        """三种合成布局总数72、基因型各36，混合笼数为0、6、12。"""
        for mix in [0., .5, 1.]:
            u = synthetic(6, mix, .5, 3252026)
            self.assertEqual(len(u), 72)
            self.assertEqual(int(u.genotype.sum()), 36)
            self.assertEqual(int(u.groupby('cage').genotype.nunique().eq(2).sum()), round(12 * mix))
            _, _, info = matrices(u)
            self.assertEqual(info['fixed_estimable'], mix > 0)
            if mix == 1:
                self.assertAlmostEqual(info['information_fraction'], 1.)

    def test_sex_in_actual_design(self):
        """真实布局保留年龄和性别，不因它们的模拟系数为零就删去。"""
        u = synthetic(6, .5, .5, 7)
        u['age'] = np.tile([4, 8, 12], 24)
        u['sex'] = np.tile(['F', 'M'], 36)
        x, _, _ = matrices(u)
        self.assertEqual(list(x), ['Intercept', 'genotype', 'age_centered', 'male'])

    def test_simulation_registry(self):
        """默认情境固定32个，重复数和种子都写入清单。"""
        u = synthetic(6, .5, .5, 7)
        u['genotype'] = u.genotype.map({0: 'WT', 1: '5X'})
        u['age'] = np.tile([4, 8, 12], 24)
        u['sex'] = np.tile(['F', 'M'], 36)
        u['SAMPLE'] = [str(200+i) for i in range(len(u))]
        (self.work / 'prepared').mkdir()
        (self.work / 'results').mkdir()
        u.to_csv(self.work / 'prepared/mouse_units.csv', index=False)
        prepare(self.work, self.config)
        tasks = pd.read_csv(self.work / 'results/simulation_scenarios.csv')
        self.assertEqual(len(tasks), 32)
        self.assertEqual(tasks.layout.nunique(), 4)
        self.assertTrue(tasks.seed.is_unique)
        self.assertTrue(tasks.planned.eq(500).all())

    def test_simulation_row_order(self):
        """同一动物清单重排行后，随机数分配使用的规范行序保持不变。"""
        u = pd.DataFrame(dict(SAMPLE=['222','211','221'], cage=['2.2','2.1','2.2'],
                              genotype=[1,0,0], age=[12,12,8], sex=['F','M','F']))
        a = canonical_actual(u)
        b = canonical_actual(u.sample(frac=1, random_state=7))
        pd.testing.assert_frame_equal(a,b)
        self.assertEqual(a.cage.tolist(),[0,1,1])

    def test_simulation_requires_sample_id(self):
        """没有来源编号不能按当前行序继续模拟，否则同一种子不再可复现。"""
        with self.assertRaisesRegex(ValueError,'SAMPLE'):
            canonical_actual(pd.DataFrame(dict(cage=['a','b'],genotype=[0,1])))

    def test_simulation_numeric_cage(self):
        """合成笼位保持整数0到11，不依赖R字符排序。"""
        u = synthetic(6,.5,.5,3252026)
        self.assertTrue(pd.api.types.is_integer_dtype(u.cage))
        self.assertEqual(sorted(u.cage.unique()),list(range(12)))

    def test_simulation_missing_sample_id(self):
        """缺失编号不能被字符串转换变成看似有效的nan标签。"""
        for value in [None, np.nan, '  ']:
            with self.assertRaisesRegex(ValueError, 'SAMPLE'):
                canonical_actual(pd.DataFrame(dict(SAMPLE=['a', value], cage=['x', 'y'])))

    def test_bundle_contains_failures_and_code(self):
        """审阅包保留小型失败明细与执行代码，省略大规模逐次模拟表。"""
        (self.work / 'results').mkdir()
        (self.work / 'results/simulation_failures.csv').write_text('reason\nerror\n')
        (self.work / 'results/simulation_replicates.csv.gz').write_bytes(b'large')
        run.bundle(self.work)
        with zipfile.ZipFile(self.work / 'review_bundle.zip') as z:
            names = set(z.namelist())
        self.assertIn('results/simulation_failures.csv', names)
        self.assertIn('executed_code/analysis/models.R', names)
        self.assertIn('executed_code/sha256.json', names)
        self.assertNotIn('results/simulation_replicates.csv.gz', names)


if __name__ == '__main__':
    unittest.main(verbosity=2)
