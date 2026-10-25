"""笼位分类标签与样本身份表的回归测试。"""
from pathlib import Path
import tempfile
import unittest
import zipfile
import numpy as np
import pandas as pd
import openpyxl
from pipeline.sources import cage_id, scalar_id, mouse_manifest, MOUSE_TABLES


class CageTests(unittest.TestCase):
    """区分样本编号限制与笼位标签，贯通GEO表到发表补充表的连接。"""

    def test_fractional_cage(self):
        """小数笼号不能截断，整数消除Excel存储带来的小数尾部。"""
        self.assertEqual(cage_id(6.5), '6.5')
        self.assertEqual(cage_id(6.0), '6')
        self.assertNotEqual(cage_id(6.5), cage_id(6))
        self.assertEqual(cage_id(' cage-A '), 'cage-A')
        self.assertEqual(cage_id('006'), '006')

    def test_invalid_cage(self):
        """缺失值、无穷值和逻辑值不得变成有效笼位。"""
        for value in (None, np.nan, np.inf, True, '', 'NA'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                cage_id(value)

    def test_sample_still_integer(self):
        """笼位可以有小数标签，SAMPLE编号仍须是整数。"""
        with self.assertRaises(ValueError):
            scalar_id(6.5)

    def test_manifest_roundtrip(self):
        """六份样本表连接S2后保留6和6.5两个笼，导出的身份一致。"""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / 'raw'; raw.mkdir()
            prepared = root / 'prepared'; prepared.mkdir()
            book = openpyxl.Workbook()
            sheet = book.active; sheet.title = 'Supp. Table 2'
            sheet.append(['SAMPLE', 'ID', 'Cage.No', 'EarNotch.No', 'Sex', 'Age', 'Genotype', 'Slide.1'])
            for i, name in enumerate(MOUSE_TABLES):
                sample = 100 + i
                pd.DataFrame([dict(SAMPLE=sample, Sex='M', Age='8mo', Condition='5X')]).to_csv(raw/name, index=False)
                sheet.append([sample, f'animal-{i}', [6, 6.5, 7, 8, 9, 10][i], 1, 'M', '8mo', '5X', 'slide-A'])
            workbook = root / 'tables.xlsx'; book.save(workbook); book.close()
            with zipfile.ZipFile(raw/'mouse_tables.zip', 'w') as archive:
                archive.write(workbook, 'tables.xlsx')
            mouse_manifest(raw, prepared)
            result = pd.read_csv(prepared/'mouse_manifest.csv', dtype={'cage': str})
            self.assertEqual(result.cage.tolist(), ['6', '6.5', '7', '8', '9', '10'])
            self.assertEqual(result.animal_id.tolist(), [f'animal-{i}' for i in range(6)])
