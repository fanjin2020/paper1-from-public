"""Numeric reporting contracts; fixtures are synthetic, never manuscript evidence."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import stats

from pipeline import reporting as R


class ReportingTests(unittest.TestCase):
    def test_geometry_and_vif(self):
        rng=np.random.default_rng(83)
        x=pd.DataFrame({'intercept':np.ones(80),'a':rng.normal(size=80),'b':rng.normal(size=80)})
        row,vif=R.geometry(x)
        self.assertEqual(row['residual_df'],77)
        self.assertTrue(all(v['VIF']>=1 for v in vif))
        x['b']=x.a
        self.assertEqual(R.geometry(x)[0]['rank'],2)

    def test_scope_missing_tests_remain_in_family(self):
        x=pd.DataFrame(dict(model=['base']*3,contrast=['clinical']*3,context=['A']*3,
                            spec=['A::clinical']*3,p=[.01,np.nan,.9]))
        z=R.scope_values(x)
        self.assertAlmostEqual(z.q_all_planned.iloc[0],.03)
        self.assertAlmostEqual(z.q_all_available.iloc[0],.02)
        self.assertTrue(pd.isna(z.q_all_planned.iloc[1]))

    @staticmethod
    def metadata():
        n=24
        return pd.DataFrame(dict(donor_id=[f'D{i}' for i in range(n)],n_cells=100,
            age=np.arange(n)+60,sex=np.tile(['male','female'],n//2),pathology=np.arange(n)%6,
            PMI=np.arange(n)/2+3,RIN=np.arange(n)/10+6,pH=np.arange(n)/100+6,
            clinical_group=['No_dementia']*(n//2)+['Dementia']*(n//2)))

    def test_donors_not_pseudoreplicated_and_bootstrap_deterministic(self):
        d=self.metadata();d=pd.concat([d,d],ignore_index=True)
        a=R.donor_tables(d,repetitions=20);b=R.donor_tables(d,repetitions=20)
        self.assertEqual(len(a[0]),24)
        self.assertEqual(len(a[3]),6)
        for x,y in zip(a,b):pd.testing.assert_frame_equal(x,y)
        d.loc[24,'age']=100
        with self.assertRaisesRegex(ValueError,'Inconsistent'):
            R.donor_tables(d,repetitions=2)

    def test_export_count_handles_quoted_newlines(self):
        with tempfile.TemporaryDirectory() as d:
            ex=R.Exports(Path(d))
            ex.add('S01','example',frame=pd.DataFrame({'a':['line1\nline2','b']}))
            self.assertEqual(ex.index[0]['rows'],2)
            with self.assertRaisesRegex(ValueError,'Duplicate'):
                ex.add('S01','example',frame=pd.DataFrame({'a':[1]}))

    def test_large_fields_roundtrip_and_restore_csv_limit(self):
        long_value = 'quoted,"value"\n' * 20000
        frame = pd.DataFrame({'setting':['jobs','seed'], 'value':[long_value,'17']})
        previous = csv.field_size_limit(131072)
        try:
            with tempfile.TemporaryDirectory() as d:
                root = Path(d)
                ex = R.Exports(root / 'export')
                plain = ex.add('S17', 'run_settings', frame=frame)
                self.assertEqual(ex.index[-1]['rows'], 2)
                self.assertEqual(pd.read_csv(plain).value.iloc[0], long_value)
                self.assertEqual(csv.field_size_limit(), 131072)
                zipped = root / 'settings.csv.gz'
                frame.to_csv(zipped, index=False)
                target = ex.add('S20', 'compressed_settings', source=zipped)
                self.assertEqual(ex.index[-1]['rows'], 2)
                self.assertEqual(target.read_bytes(), zipped.read_bytes())
                self.assertEqual(csv.field_size_limit(), 131072)
                with self.assertRaises(FileNotFoundError):
                    R.row_count(root / 'missing.csv')
                self.assertEqual(csv.field_size_limit(), 131072)
        finally:
            csv.field_size_limit(previous)

    def test_block_map_rejects_duplicate_semantic_owners(self):
        with tempfile.TemporaryDirectory() as d:
            ex=R.Exports(Path(d))
            expected=[f'S{s}.{i}' for s,n in enumerate([3,3,11,4,5,5,4,3,4,1,6,11,4,5,4,4,7,3,6,12],1) for i in range(1,n+1)]
            ex.index=[dict(section='S01',blocks=','.join(expected),file='all.csv',description='test')]
            self.assertEqual(len(ex.finish()),1)
            ex.index.append(dict(section='S01',blocks='S1.1',file='second.csv',description='ambiguous'))
            with self.assertRaisesRegex(ValueError,'explicit combined'):
                ex.finish()


if __name__=='__main__':
    unittest.main()
