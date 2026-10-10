"""结果完整性审计测试；表内数值为人工构造，不是研究输出。"""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from pipeline.audit import assess
from pipeline.io import dump


class AuditTests(unittest.TestCase):
    """用最小的完整结果网格验证失败分类和核对规则。"""

    def setUp(self):
        """保留实际网格数量，但每个规格只放一个人工基因。"""
        self.temp = tempfile.TemporaryDirectory()
        self.work = Path(self.temp.name)
        self.r = self.work / 'results'; self.r.mkdir()
        self.specs = ['whole_linear', 'whole_age_factor', 'whole_age_factor_slide',
                      'cortex_linear', 'cortex_age_factor', 'cortex_age_factor_slide']
        self.models = ['baseline', 'cage_fixed', 'cage_correlation', 'dream_shared', 'dream_native']
        dump(self.work / 'run_config.json', dict(simulation_repeats=500))
        pd.DataFrame({'gene_symbol': ['A']}).to_csv(self.r / 'panel.csv', index=False)
        pd.DataFrame({'gene_symbol': ['A']}).to_csv(self.r / 'mouse_background.csv', index=False)
        h = pd.DataFrame([dict(gene_symbol='A', spec='Ast::clinical_status', context='Ast',
             contrast='clinical_status', model=m, beta=b, SE=.1, p=p, q_joint=p)
             for m, b, p in [('base', 1., .01), ('quality', .9, .2)]])
        h.to_csv(self.r / 'human_results.csv.gz', index=False)
        mouse = pd.DataFrame([dict(gene_symbol='A', spec=s, model=m) for s in self.specs for m in self.models])
        mouse.to_csv(self.r / 'mouse_results.csv.gz', index=False)
        status = mouse.assign(planned=1, available=1, failed=0)
        status.to_csv(self.r / 'mouse_model_status.csv', index=False)
        pd.DataFrame(columns=['gene_symbol','spec','model','status']).to_csv(self.r / 'mouse_fit_failures.csv', index=False)
        var = pd.DataFrame([dict(gene_symbol='A', spec=s, weight_policy=w, status='available')
                            for s in self.specs for w in ['shared', 'native']])
        var.to_csv(self.r / 'mouse_variance.csv.gz', index=False)
        prec = pd.DataFrame([dict(gene_symbol='A', spec=s, SE_ratio=1.5, design_factor=1.5,
                                 posterior_factor=1., reweight_factor=1.) for s in self.specs])
        prec.to_csv(self.r / 'mouse_precision.csv.gz', index=False)
        impl = pd.DataFrame([dict(gene_symbol='A', spec=s, branch=b, unscaled_exact=.1, unscaled_direct=.1,
                                 unscaled_contrast_exact=.1,unscaled_contrast_refit=.1)
                             for s in self.specs for b in ['baseline', 'cage_fixed']])
        impl.to_csv(self.r / 'mouse_SE_implementation_audit.csv.gz', index=False)
        pd.DataFrame(columns=['gene_symbol','status','message','spec','model']).to_csv(
            self.r/'mouse_dream_diagnostics.csv',index=False)
        sim = pd.DataFrame([dict(scenario=f's{i:03d}', model=m, planned=500, available=500,
                                numerical_failed=0, not_estimable=0)
                            for i in range(32) for m in self.models[:3]+['dream']])
        sim.to_csv(self.r / 'simulation_summary.csv', index=False)
        pd.DataFrame(columns=['scenario','model','status','reason']).to_csv(self.r / 'simulation_failures.csv', index=False)

    def tearDown(self):
        """清理本测试的临时表，不触及任何实际结果目录。"""
        self.temp.cleanup()

    def test_complete_but_different_numbers(self):
        """旧稿数字不同不修改估计，只单独记录差异。"""
        result = assess(self.work)
        self.assertEqual(result['status'], 'completed')
        self.assertGreater(result['manuscript_anchor_differences'], 0)
        self.assertFalse(result['historical_manuscript_reproduction_certified'])

    def test_numerical_failure_changes_status(self):
        """一次数值失败必须影响总状态，不能被其他成功结果抵消。"""
        p = self.r / 'simulation_summary.csv'
        s = pd.read_csv(p); s.loc[0, ['available', 'numerical_failed']] = [499, 1]
        s.to_csv(p, index=False)
        pd.DataFrame([dict(scenario='s000', model='baseline', status='failed', reason='test')]).to_csv(
            self.r / 'simulation_failures.csv', index=False)
        result = assess(self.work)
        self.assertEqual(result['status'], 'completed_with_fit_failures')
        self.assertEqual(result['failures']['simulation_fits'], 1)

    def test_nonestimability_is_separate(self):
        """结构性不可估计不冒充数值失败。"""
        p = self.r / 'simulation_summary.csv'
        s = pd.read_csv(p); s.loc[1, ['available', 'not_estimable']] = [0, 500]
        s.to_csv(p, index=False)
        pd.DataFrame([dict(scenario='s000', model='cage_fixed', status='not_estimable', reason='test')]*500).to_csv(
            self.r / 'simulation_failures.csv', index=False)
        result = assess(self.work)
        self.assertEqual(result['status'], 'completed')
        self.assertEqual(result['simulation_not_estimable'], 500)

    def test_failure_denominator(self):
        """计划数与可用、失败、不可估计之和不同就必须停止。"""
        p = self.r / 'simulation_summary.csv'
        s = pd.read_csv(p); s.loc[0, 'available'] = 499; s.to_csv(p, index=False)
        with self.assertRaisesRegex(ValueError, '分母不能闭合'):
            assess(self.work)

    def test_missing_native_variance(self):
        """只输出共享权重ICC时不能通过完整性检查。"""
        p = self.r / 'mouse_variance.csv.gz'
        s = pd.read_csv(p); s.loc[s.weight_policy.eq('shared')].to_csv(p, index=False)
        with self.assertRaisesRegex(ValueError, '方差组分缺少'):
            assess(self.work)

    def test_inexact_standard_error(self):
        """直接SE与精确协方差不一致时，不允许结束为成功。"""
        p = self.r / 'mouse_SE_implementation_audit.csv.gz'
        s = pd.read_csv(p); s.loc[0, 'unscaled_direct'] = .2; s.to_csv(p, index=False)
        with self.assertRaisesRegex(ValueError, '精确协方差'):
            assess(self.work)

    def test_unresolved_warning(self):
        """有限估计不能掩盖没有定位的收敛警告。"""
        pd.DataFrame([dict(gene_symbol=None,status='warning_not_reproduced',message='negative eigenvalue')]).to_csv(
            self.r/'mouse_dream_diagnostics.csv',index=False)
        result=assess(self.work)
        self.assertEqual(result['status'],'completed_with_unresolved_warnings')
        self.assertEqual(result['dream_unresolved_warnings'],1)

    def test_exact_contrast_check(self):
        """原参数化精确对比不等于直接系数时，不能只凭两者系数相同放行。"""
        p=self.r/'mouse_SE_implementation_audit.csv.gz'
        s=pd.read_csv(p);s.loc[0,'unscaled_contrast_exact']=.3;s.to_csv(p,index=False)
        with self.assertRaisesRegex(ValueError,'精确协方差'):
            assess(self.work)


if __name__ == '__main__':
    unittest.main(verbosity=2)
