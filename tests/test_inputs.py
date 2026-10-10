"""下载缓存、标识符、矩阵聚合及恢复机制的回归测试。"""
import ast
import gzip
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from pipeline.io import download, dump, sha
from pipeline.sources import scalar_id, age_month, genotype, read_geo_matrix
from pipeline.seaad import collapse_genes, combine, count_sum_sparse, gene_symbols


class InputTests(unittest.TestCase):
    """用独立的小型数据触发正常分支和拒绝分支，不下载真实大对象。"""

    def test_ids(self):
        """编号仅删除整数的小数尾部，不把动物标签中的数字猜作样本号。"""
        self.assertEqual(scalar_id(342.0), '342')
        with self.assertRaises(ValueError):
            scalar_id('animal342')

    def test_age(self):
        """月份单位可解析，日期字符串不作为年龄。"""
        self.assertEqual(age_month('8mo'), 8)
        self.assertTrue(np.isnan(age_month('12/10/20')))

    def test_genotype(self):
        """条件缩写只按明确字典映射。"""
        self.assertEqual(genotype('5xFAD'), '5X')
        self.assertIsNone(genotype('AD'))

    def test_gene_sum(self):
        """同名计数行相加，不使用均值或只保留首行。"""
        names = np.array(['GFAP','GFAP']+[f'GENE{i}' for i in range(101)])
        genes, values = collapse_genes(names,np.ones((103,2)))
        np.testing.assert_array_equal(values[genes=='GFAP'],[[2,2]])

    def test_human_symbol_case(self):
        """人类orf符号按统一大小写匹配，不因来源写法漏掉同一个端点。"""
        raw = ['C15orf41', 'C16orf96', 'C1orf112', 'C3orf79', 'C4orf48']
        np.testing.assert_array_equal(gene_symbols(raw), [x.upper() for x in raw])
        self.assertEqual(gene_symbols(['ENSG000001'])[0], '')

    def test_matrix_identity(self):
        """计数列与元数据顺序不一致时必须停止。"""
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaises(ValueError):
                combine([(np.array(['A']),np.array([[1]]),['u1'])],
                        [dict(unit_id='u2')],Path(folder),'test')

    def test_csr_grouping(self):
        """稀疏分组求和与手算一致，未选中的细胞不进入聚合。"""
        group={'indptr':np.array([0,2,3,5]),'indices':np.array([0,2,1,0,1]),'data':np.array([1,2,4,5,6])}
        x=count_sum_sparse(group,[np.array([0,2])],['Ast'],3,'rows')
        np.testing.assert_array_equal(x,np.array([[6],[6],[2]]))

    def test_csr_fractional(self):
        """标准化后的非整数层不能作为UMI计数进入分析。"""
        group={'indptr':np.array([0,1]),'indices':np.array([0]),'data':np.array([1.2])}
        with self.assertRaises(ValueError):
            count_sum_sparse(group,[np.array([0])],['Ast'],1,'rows')

    def test_csr_duplicate_cells(self):
        """一个核被分给两个类群时停止，不能重复贡献计数。"""
        group={'indptr':np.array([0,1]),'indices':np.array([0]),'data':np.array([1])}
        with self.assertRaises(ValueError):
            count_sum_sparse(group,[np.array([0]),np.array([0])],['Ast','Oli'],1,'rows')

    def test_unregistered_cache(self):
        """复制来的旧文件没有下载收据，不能成为新分析输入。"""
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'counts.gz';p.write_bytes(b'old')
            with self.assertRaises(ValueError):
                download('https://example.invalid/a',p)

    def test_changed_cache(self):
        """缓存内容被修改后，即使长度相同也不能复用。"""
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'file';p.write_bytes(b'first')
            dump(p.with_name('file.download.json'),dict(url='test',bytes=5,sha256=sha(p)))
            p.write_bytes(b'other')
            with self.assertRaises(ValueError):
                download('test',p)

    def test_partial_http(self):
        """分段响应尚未达到Content-Range总长时必须继续，不登记为完整下载。"""
        class Response(io.BytesIO):
            """替代HTTP响应，模拟服务端分两段返回同一文件。"""
            def __init__(self, data, content_range):
                """设置续传校验所需的响应头。"""
                super().__init__(data)
                self.status=206
                self.headers={'Content-Length':'3','Content-Range':content_range,'ETag':'same'}
        responses=[Response(b'abc','bytes 0-2/6'),Response(b'def','bytes 3-5/6')]
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'file'
            with patch('urllib.request.urlopen',side_effect=responses),patch('time.sleep'):
                download('https://example.test/file',p)
            self.assertEqual(p.read_bytes(),b'abcdef')

    def test_geo_parser(self):
        """GEO表达列和样本特征按同一编号连接。"""
        text='!Sample_geo_accession\t"GSM1"\t"GSM2"\n!Sample_characteristics_ch1\t"individual: A"\t"individual: B"\n!series_matrix_table_begin\n"ID_REF"\t"GSM1"\t"GSM2"\n"p1"\t1\t2\n!series_matrix_table_end\n'
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'matrix.gz'
            with gzip.open(p,'wt') as f:
                f.write(text)
            x,m=read_geo_matrix(p)
            self.assertEqual(x.loc['p1','GSM2'],2)
            self.assertEqual(m.individual.tolist(),['A','B'])

    def test_truncated_gzip(self):
        """下载不完整的gzip即使前部能读，也不能通过整个矩阵解析。"""
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'broken.gz';p.write_bytes(gzip.compress(b'abc')[:-5])
            with self.assertRaises((EOFError,ValueError)):
                read_geo_matrix(p)

    def test_sources_are_public(self):
        """发布包不包含用户私有绝对路径或历史分析目录引用。"""
        root=Path(__file__).resolve().parents[1]
        for folder in ('pipeline','analysis'):
            for p in (root/folder).glob('*'):
                if p.suffix in ('.py','.R'):
                    text=p.read_text(encoding='utf-8')
                    self.assertNotIn('/home/research/',text)
                    self.assertNotIn('paper1_public_evidence',text)

    def test_docstrings(self):
        """所有Python函数和类保留说明，便于读者审阅和维护。"""
        root=Path(__file__).resolve().parents[1]
        for p in [root/'run.py',root/'run_dream_repair.py',*(root/'pipeline').glob('*.py')]:
            for node in ast.walk(ast.parse(p.read_text(encoding='utf-8'))):
                if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef)):
                    self.assertTrue(ast.get_docstring(node),f'{p.name}:{node.lineno}缺少说明')


if __name__=='__main__':
    unittest.main(verbosity=2)
