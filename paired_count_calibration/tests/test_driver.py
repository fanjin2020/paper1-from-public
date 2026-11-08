import copy
import json
from pathlib import Path
import tempfile
import unittest
from paired_count_calibration import driver as d


def fixture(out):
    p=dict(output=str(out),run_id='test',settings=dict(genes=200,reps=2,seed=1024),scenarios=d.scenarios('smoke'))
    s=p['scenarios'][0]
    rows=[]
    for mode,test in sorted(d.KEYS):
        r=dict(weight_mode=mode,test=test,genes=200,retained=180,available=170,filtered=20,
               unavailable_retained=10,paired_truth=0,p_histogram=[170]+[0]*19,
               coverage95=.95,delivered_coverage95=.8,fpr_010=.01,fpr_050=.05,fpr_100=.1)
        for family in ('planned','retained'):
            for level in ('010','050','100'):
                r[f'discoveries_{family}_{level}']=0
                r[f'fdp_{family}_{level}']=0
        rows.append(r)
    result=dict(run_id='test',scenario=s['id'],replicate=1,seed=1024+s['index']*100000+1,
                metrics=rows,checks={m:dict(lmfit_max_difference=1e-14) for m in ('shared','design_specific')})
    return p,s,result


class DriverTests(unittest.TestCase):
    def test_full_grid(self):
        a=d.scenarios('full')
        self.assertEqual(len(a),24)
        self.assertEqual(len({x['id'] for x in a}),24)
        self.assertEqual(len({(x['n'],x['rho'],x['depth'],x['exposure_fraction']) for x in a}),24)

    def test_smoke_grid(self):
        a=d.scenarios('smoke')
        self.assertEqual({x['depth'] for x in a},{'equal','iid_lognormal'})
        self.assertEqual({x['rho'] for x in a},{.3,.7})
        self.assertEqual({x['exposure_fraction'] for x in a},{0,.2})

    def test_no_duplicate_seeds(self):
        seeds={20261107+s['index']*100000+r for s in d.scenarios('full') for r in range(1,1001)}
        self.assertEqual(len(seeds),24000)

    def test_receipt_required_not_job(self):
        with tempfile.TemporaryDirectory() as t:
            p,s,r=fixture(t);f=d.folder(p,s,1)
            d.atomic_json(f/'job.json',{})
            self.assertFalse(d.valid_receipt(p,s,1))
            d.atomic_json(f/'result.json',r)
            self.assertFalse(d.valid_receipt(p,s,1))
            d.publish_receipt(p,s,1)
            self.assertTrue(d.valid_receipt(p,s,1))
            self.assertEqual(d.progress(p)[0],1)

    def test_corruption_is_pending(self):
        with tempfile.TemporaryDirectory() as t:
            p,s,r=fixture(t);f=d.folder(p,s,1)
            d.atomic_json(f/'result.json',r);d.publish_receipt(p,s,1)
            r['metrics'][0]['fpr_050']=.9
            d.atomic_json(f/'result.json',r)
            self.assertFalse(d.valid_receipt(p,s,1))
            self.assertEqual(d.progress(p)[0],0)

    def test_invalid_metrics(self):
        p,s,r=fixture('not-used')
        cases=[]
        a=copy.deepcopy(r);a['metrics'].pop();cases.append(a)
        a=copy.deepcopy(r);a['metrics'][0]['available']=201;cases.append(a)
        a=copy.deepcopy(r);a['metrics'][0]['fpr_050']=None;cases.append(a)
        a=copy.deepcopy(r);a['metrics'][0]['p_histogram']=[0]*20;cases.append(a)
        a=copy.deepcopy(r);a['metrics'][0]['fdp_planned_050']=1;cases.append(a)
        a=copy.deepcopy(r);a['checks']['shared']['lmfit_max_difference']=.1;cases.append(a)
        for a in cases:
            with self.assertRaises(ValueError):d.validate_result(a,p,s,1)

    def test_bad_seed(self):
        p,s,r=fixture('not-used');r['seed']+=1
        with self.assertRaises(ValueError):d.validate_result(r,p,s,1)

    def test_source_identity(self):
        self.assertEqual(d.identity(dict(a=1,b=2)),d.identity(dict(b=2,a=1)))
        self.assertNotEqual(d.identity(dict(a=1)),d.identity(dict(a=2)))

    def test_mcse_dataset_unit(self):
        avg,se=d.mean_mc([0,1,0,1])
        self.assertEqual(avg,.5)
        self.assertAlmostEqual(se,1/(12**.5))
        self.assertIsNone(d.mean_mc([1])[1])

    def test_wilson_boundaries(self):
        lo,hi=d.wilson(0,1000)
        self.assertEqual(lo,0)
        self.assertGreater(hi,0)
        self.assertLess(hi,.004)
        lo,hi=d.wilson(50,1000)
        self.assertLess(lo,.05);self.assertGreater(hi,.05)

    def test_atomic_json_rejects_nan(self):
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaises(ValueError):d.atomic_json(Path(t)/'test.json',dict(bad=float('nan')))

    def test_writer_lock(self):
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)
            with d.output_lock(out):self.assertTrue(d.writer_active(out))
            self.assertFalse(d.writer_active(out))


if __name__=='__main__':unittest.main()
