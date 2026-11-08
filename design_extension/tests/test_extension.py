import copy
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from design_extension.geometry import design_diagnostics, linear_operator, residual
from design_extension.study import create_designs, scenarios, gaussian_job, bh
from design_extension.report import fit_rule, choose, summarize, csv_rows, paired_mean_shift_contrasts
from design_extension.driver import sha, identity, valid_receipt, atomic_json, validate_metrics, validate_count_health, check_output_location


def layout():
    g=np.repeat(np.arange(12),4)
    x=np.concatenate([np.tile([0.,1.],2) for _ in range(6)]+[np.zeros(4)]*3+[np.ones(4)]*3)
    return x,g


class GeometryTests(unittest.TestCase):
    def test_within_ratio_and_r_squared(self):
        x,g=layout();z=np.eye(12)[g]
        d=design_diagnostics(x,g)
        self.assertAlmostEqual(d['ewi'],.5)
        self.assertAlmostEqual(d['conditional_ewi'],.5)
        r2=1-float(residual(x,z)@residual(x,z))/np.sum((x-x.mean())**2)
        self.assertAlmostEqual(d['ewi'],1-r2)

    def test_pure_balanced_and_covariate_collinear(self):
        g=np.repeat(np.arange(12),4)
        x=np.repeat([0.,1.]*6,4)
        d=design_diagnostics(x,g)
        self.assertFalse(d['within_estimable'])
        self.assertAlmostEqual(d['ewi'],0.)
        self.assertIsNone(linear_operator(x,g,model='fixed'))
        self.assertIsNone(linear_operator(x,g,model='mundlak_known_covariance'))
        x=np.tile([0.,1.],24)
        self.assertAlmostEqual(design_diagnostics(x,g)['ewi'],1.)
        self.assertFalse(design_diagnostics(x,g,x[:,None])['baseline_estimable'])

    def test_analytic_bias_and_group_bound(self):
        x,g=layout();z=np.eye(12)[g];means=z@(z.T@x/z.sum(0))-x.mean()
        op=linear_operator(x,g)
        self.assertAlmostEqual(op['gamma_bias_multiplier'],1-design_diagnostics(x,g)['ewi'])
        rng=np.random.default_rng(42);y=rng.normal(size=len(x));shift=rng.normal(size=12)
        for m in ('baseline','fixed','re_known_covariance','mundlak_known_covariance'):
            op=linear_operator(x,g,model=m,cage_sd=.3)
            delta=op['a']@(y+.7*means)-op['a']@y
            self.assertAlmostEqual(delta,.7*op['gamma_bias_multiplier'])
            self.assertLessEqual(abs(op['a']@z@shift),op['group_shift_l2_multiplier']*np.linalg.norm(shift)+1e-12)
        self.assertAlmostEqual(linear_operator(x,g,model='fixed')['gamma_bias_multiplier'],0.)

    def test_mundlak_coefficient_identity(self):
        x,g=layout();rng=np.random.default_rng(331)
        cov=rng.normal(size=(len(x),2));y=rng.normal(size=(len(x),8))
        fixed=linear_operator(x,g,cov,'fixed')['a']@y
        for sd in (0.,.3,2.):
            cre=linear_operator(x,g,cov,'mundlak_known_covariance',sd)['a']@y
            np.testing.assert_allclose(fixed,cre,atol=1e-11)

    def test_row_group_and_column_scaling_invariance(self):
        x,g=layout();rng=np.random.default_rng(911);cov=rng.normal(size=(len(x),2))
        d=design_diagnostics(x,g,cov);order=rng.permutation(len(x))
        e=design_diagnostics(x[order],np.array(['c'+str(v) for v in g[order]]),cov[order]*[1e-5,1e5])
        for key in ('ewi','conditional_ewi','within_information','effective_information_groups'):
            self.assertAlmostEqual(d[key],e[key],places=9)

    def test_invalid_design(self):
        with self.assertRaises(ValueError):design_diagnostics([0.,1.],['a',None])
        with self.assertRaises(ValueError):design_diagnostics([0.,np.nan],['a','b'])
        with self.assertRaises(ValueError):design_diagnostics([0.],['a'])


class SimulationAndRuleTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.designs=create_designs(self.root/'designs',321)

    def tearDown(self):self.tmp.cleanup()

    def test_grid_and_holdout(self):
        self.assertEqual(len(self.designs),12)
        self.assertEqual(sum(d['split']=='holdout' for d in self.designs),4)
        self.assertTrue(all(d['groups']==18 for d in self.designs if d['split']=='holdout'))
        self.assertEqual(len(scenarios('full')),8)

    def test_gaussian_validation_and_repeatability(self):
        d=self.designs[1]
        job=dict(id='test',phase='gaussian',design=d['design'],design_path=d['path'],split=d['split'],
                 design_index=1,start=0,end=4,genes=100,scenario=scenarios('pilot')[1],seed=321)
        gaussian_job(job,self.root/'a');gaussian_job(job,self.root/'b')
        a=csv_rows(self.root/'a/metrics.csv');b=csv_rows(self.root/'b/metrics.csv')
        self.assertEqual(a,b);validate_metrics(a,job)
        damaged=copy.deepcopy(a);damaged[0]['scenario']='9'
        with self.assertRaises(ValueError):validate_metrics(damaged,job)
        damaged=copy.deepcopy(a);damaged[0]['FDP']='1.2'
        with self.assertRaises(ValueError):validate_metrics(damaged,job)

    def test_zero_discoveries_retain_uncertainty(self):
        _,upper,n=summarize([{'FDP':0}]*50,'FDP')
        self.assertEqual(n,50);self.assertGreater(upper,0.05)
        self.assertLess(upper,.075)

    def test_bh_planned_family(self):
        np.testing.assert_allclose(bh([.01,.04,np.nan])[:2],[.03,.06])
        self.assertTrue(np.isnan(bh([.01,np.nan])[1]))

    def test_paired_mean_shift_summary(self):
        s=scenarios('pilot');rows=[]
        for scenario in s:
            for rep in range(3):
                rows.append(dict(design='test',split='training',scenario=scenario['scenario'],model='baseline',
                                 replicate=rep,FDP=.1*rep+scenario['cage_mean_amplitude']))
        result=paired_mean_shift_contrasts(rows,s)
        self.assertEqual(len(result),2)
        for row in result:
            self.assertAlmostEqual(row['mean_FDP_difference'],.5)
            self.assertAlmostEqual(row['paired_mcse'],0.)

    def test_entire_model_failure_is_fatal(self):
        d=self.designs[1]
        rows=[dict(model=m,availability=0. if m=='dream_native' else 1.) for m in
              ('baseline','cage_fixed','cage_correlation','dream_native','mundlak_native')]
        job=dict(design_path=d['path'],output='test')
        with self.assertRaises(ValueError):validate_count_health(rows,job)
        rows[3]['availability']=1.;validate_count_health(rows,job)

    def test_pure_design_failures_are_structural(self):
        d=self.designs[0]
        rows=[dict(model=m,availability=0. if m in ('cage_fixed','mundlak_native') else 1.) for m in
              ('baseline','cage_fixed','cage_correlation','dream_native','mundlak_native')]
        validate_count_health(rows,dict(design_path=d['path'],output='test'))

    def training_rows(self):
        return [dict(phase='counts',design=d['design'],split=d['split'],scenario=0,model='cage_fixed',
                 usable_datasets=50,availability=1.,fdr_one_sided_95_mc_upper=.06,power=.6,empirical_fdr=.04)
                for d in self.designs]

    def test_holdout_results_never_fit_rule(self):
        rows=self.training_rows();before=fit_rule(rows,self.designs)
        for row in rows:
            if row['split']=='holdout':row['power']=999;row['empirical_fdr']=999
        self.assertEqual(before,fit_rule(rows,self.designs))
        d=next(d for d in self.designs if d['split']=='holdout' and d['within_estimable'])
        self.assertEqual(choose(before,d,0)['model'],'cage_fixed')

    def test_abstain_pure_and_pilot(self):
        rows=self.training_rows();pure=next(d for d in self.designs if not d['within_estimable'])
        self.assertIsNone(choose(fit_rule(rows,self.designs),pure,0)['model'])
        rule=fit_rule(rows,self.designs,validation_only=True)
        self.assertIsNone(choose(rule,self.designs[1],0)['model'])
        for row in rows:row['usable_datasets']=3
        for d in self.designs:self.assertIsNone(choose(fit_rule(rows,self.designs),d,0)['model'])


class ReceiptTests(unittest.TestCase):
    def test_output_never_overlaps_source_or_inputs(self):
        with tempfile.TemporaryDirectory() as root:
            root=Path(root);code=root/'code';data=root/'data'
            for forbidden in (code,code/'run',root,data,data/'run'):
                with self.assertRaises(ValueError):check_output_location(forbidden,code,data)
            check_output_location(root/'output',code,data)

    def test_content_corruption_and_incomplete_file_map(self):
        with tempfile.TemporaryDirectory() as root:
            out=Path(root);job={'id':'j1','phase':'gaussian'};folder=out/'tasks/j1';folder.mkdir(parents=True)
            for n in ('metrics.csv','sensitivity.csv'):(folder/n).write_text('a,b\n1,2\n')
            r=dict(run_id='x',job_hash=identity(job),files={n:dict(bytes=(folder/n).stat().st_size,
                   sha256=sha(folder/n)) for n in ('metrics.csv','sensitivity.csv')})
            atomic_json(folder/'receipt.json',r);self.assertTrue(valid_receipt(out,job,'x'))
            (folder/'metrics.csv').write_text('a,b\n2,1\n')
            self.assertFalse(valid_receipt(out,job,'x'))
            del r['files']['metrics.csv'];atomic_json(folder/'receipt.json',r)
            self.assertFalse(valid_receipt(out,job,'x'))


if __name__=='__main__':unittest.main()
