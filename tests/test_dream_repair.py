"""定点修复的合并、校正分母和来源核验测试。"""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from pipeline.io import dump, sha
from run_dream_repair import bh_planned, merge_repaired, verify_source


class DreamRepairTests(unittest.TestCase):
    """模拟一个恢复基因和一个无关模型，验证修复不会丢行或改动无关估计。"""

    def test_bh_missing_denominator(self):
        """两个有效检验加一个缺失位置时，BH仍按三项校正。"""
        q = bh_planned([.01,.04,np.nan])
        np.testing.assert_allclose(q[:2],[.03,.06])
        self.assertTrue(np.isnan(q[2]))

    def test_bh_invalid(self):
        """非法p值不能通过校正被截断成看似正常的结果。"""
        with self.assertRaises(ValueError):
            bh_planned([-.1,.5])

    def test_merging(self):
        """完整替换受影响模型，非目标模型的系数不变，计划校正范围不变。"""
        with tempfile.TemporaryDirectory() as tmp:
            w = Path(tmp); r = w/'results'; r.mkdir()
            rows = [dict(gene_symbol=g, spec='whole_age_factor',model=m,beta=1.,SE=.2,
                    df=30.,p=.01,low=.6,high=1.4,status='available')
                    for m in ['baseline','dream_native'] for g in ['a','b']]
            bg = pd.DataFrame(rows)
            bg.to_csv(r/'mouse_background_fits.csv.gz',index=False)
            panel = bg.copy(); panel['mouse_gene'] = panel.gene_symbol
            panel.to_csv(r/'mouse_results.csv.gz',index=False)
            replacement = bg.loc[bg.model.eq('dream_native')].copy()
            replacement['beta'] = 2.
            replacement.to_csv(r/'dream_repair_replacement_rows.csv.gz',index=False)
            d = pd.DataFrame([dict(spec='whole_age_factor',model='dream_native',gene_symbol='a',
                                  status='recovered_optimizer_agreement',message='checked')])
            for name in ['mouse_dream_diagnostics.csv','dream_repair_diagnostics.csv',
                         'mouse_model_messages.csv','dream_repair_messages.csv']:
                d.to_csv(r/name,index=False)
            targets = pd.MultiIndex.from_tuples([('whole_age_factor','dream_native')])
            merge_repaired(w,targets)
            final = pd.read_csv(r/'mouse_background_fits.csv.gz')
            self.assertEqual(len(final),4)
            self.assertTrue(final.loc[final.model.eq('baseline'),'beta'].eq(1.).all())
            self.assertTrue(final.loc[final.model.eq('dream_native'),'beta'].eq(2.).all())
            self.assertEqual(len(pd.read_csv(r/'mouse_fit_failures.csv')),0)
            self.assertEqual(len(pd.read_csv(r/'mouse_results.csv.gz')),4)
            replacement.iloc[:1].to_csv(r/'dream_repair_replacement_rows.csv.gz',index=False)
            with self.assertRaisesRegex(ValueError,'全部背景'):
                merge_repaired(w,targets)

    def test_source_receipts(self):
        """缺少收据或计数摘要改变时停止，不能用不完整审阅包冒充原工作目录。"""
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)
            dump(source/'run_config.json',dict(analysis_version='1.2.0'))
            dump(source/'status.json',dict(status='completed_with_fit_failures',audit=dict(integrity_checks='passed')))
            names = ['results/mouse_fit_failures.csv','results/mouse_dream_diagnostics.csv',
                     'prepared/mouse_units.csv','prepared/mouse_whole_counts.csv.gz','prepared/mouse_cortex_counts.csv.gz']
            outputs = {}
            for name in names:
                p = source/name; p.parent.mkdir(parents=True,exist_ok=True)
                if name.endswith('.csv') and 'results/' in name:
                    pd.DataFrame([dict(spec='whole_age_factor',model='dream_native',gene_symbol='Micos10',
                                       status='numerical_warning')]).to_csv(p,index=False)
                else:
                    p.write_bytes(b'fixture')
                outputs[name] = sha(p)
            for stage in ['panel','orthologs','human','mouse','simulation','calibration','prepare_mouse']:
                dump(source/'stages'/f'{stage}.json',dict(stage=stage,outputs=outputs))
            _, verified, failures = verify_source(source)
            self.assertEqual(len(verified),5)
            self.assertEqual(len(failures),1)
            (source/'prepared/mouse_whole_counts.csv.gz').write_bytes(b'altered')
            with self.assertRaisesRegex(ValueError,'已改变'):
                verify_source(source)


if __name__ == '__main__':
    unittest.main()
