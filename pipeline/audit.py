"""完成后核查结果分母、拟合失败和SE恒等式，单独输出旧稿数字的比较。"""
import json
import numpy as np
import pandas as pd
from .io import dump, require


def clinical_pairs(h):
    """只取临床状态对比，防止把40条与条件临床的6条合并解释。"""
    keys = ['gene_symbol', 'spec', 'context', 'contrast']
    h = h.loc[h.contrast.eq('clinical_status')]
    return h.loc[h.model.eq('base')].merge(h.loc[h.model.eq('quality')],
            on=keys, suffixes=('_base', '_quality'), validate='one_to_one')


def assess(work):
    """审阅的是计划完整性和实际输出；旧稿计数不同只提示，不调整结果。"""
    r = work / 'results'
    h = pd.read_csv(r / 'human_results.csv.gz')
    pair = clinical_pairs(h)
    selected = pair.loc[pair.q_joint_base.lt(.05)].copy()
    mouse = pd.read_csv(r / 'mouse_results.csv.gz')
    status = pd.read_csv(r / 'mouse_model_status.csv')
    var = pd.read_csv(r / 'mouse_variance.csv.gz')
    precision = pd.read_csv(r / 'mouse_precision.csv.gz')
    impl = pd.read_csv(r / 'mouse_SE_implementation_audit.csv.gz')
    dream_diagnostics = pd.read_csv(r / 'mouse_dream_diagnostics.csv')
    sim = pd.read_csv(r / 'simulation_summary.csv')
    panel = pd.read_csv(r / 'panel.csv')
    bg = pd.read_csv(r / 'mouse_background.csv')
    require(not h.duplicated(['gene_symbol', 'spec', 'model']).any(), '人类结果有重复计划位置')
    require(len(h) == len(panel) * h.spec.nunique() * 2, '人类计划位置不完整')
    require(len(mouse) == len(panel) * 6 * 5, '小鼠计划位置不完整')
    require(not mouse.duplicated(['gene_symbol', 'spec', 'model']).any(), '小鼠结果有重复计划位置')
    require(len(status) == 30 and status.planned.eq(len(bg)).all(), '全背景模型分母不完整')
    require(len(var) == len(bg) * 6 * 2 and set(var.weight_policy) == {'shared', 'native'}, '方差组分缺少一套权重或计划位置')
    require(not var.duplicated(['gene_symbol', 'spec', 'weight_policy']).any(), '方差组分有重复位置')
    require(len(precision) == len(bg) * 6, 'SE分解背景位置不完整')
    product = precision.design_factor * precision.posterior_factor * precision.reweight_factor
    require(np.allclose(product, precision.SE_ratio, rtol=1e-8, atol=1e-10), 'SE分解恒等式失败')
    require(len(impl) == len(bg) * 6 * 2, 'SE实现诊断缺少计划位置')
    error = np.max(np.abs(impl.unscaled_direct / impl.unscaled_exact - 1))
    require(np.isfinite(error) and error < 1e-7, '直接系数SE没有通过精确协方差检查')
    for col in ['unscaled_contrast_exact', 'unscaled_contrast_refit']:
        require(np.allclose(impl[col], impl.unscaled_exact, rtol=1e-7, atol=1e-10),
                f'{col}没有通过精确协方差检查')
    require(len(sim) == 32 * 4 and sim.scenario.nunique() == 32, '四布局模拟网格不完整')
    require((sim.available + sim.numerical_failed + sim.not_estimable == sim.planned).all(), '模拟失败分母不能闭合')
    require(sim.loc[~sim.model.eq('cage_fixed'), 'not_estimable'].eq(0).all(), '非固定效应模型出现未定义的结构性不可估计标签')
    require(not sim.duplicated(['scenario', 'model']).any(), '模拟汇总有重复位置')
    with (work / 'run_config.json').open(encoding='utf-8') as f:
        config = json.load(f)
    require(sim.planned.eq(config['simulation_repeats']).all(), '模拟重复数与配置不同')
    sf = pd.read_csv(r / 'simulation_failures.csv')
    require(len(sf) == int((sim.numerical_failed + sim.not_estimable).sum()), '模拟失败明细与汇总不一致')
    mf = pd.read_csv(r / 'mouse_fit_failures.csv')
    require(len(mf) == int(status.failed.sum()), '模型失败明细与汇总不一致')
    actual = {'human_clinical_base_significant': int(pair.q_joint_base.lt(.05).sum()),
              'human_clinical_quality_significant': int(pair.q_joint_quality.lt(.05).sum()),
              'human_clinical_abs_beta_ratio': (selected.beta_quality / selected.beta_base).abs().median(),
              'human_clinical_SE_ratio': (selected.SE_quality / selected.SE_base).median()}
    # 这些数字仅是旧稿核对锚点，不参与拟合、筛选或校正。
    expected = {'human_clinical_base_significant': (40, 0), 'human_clinical_quality_significant': (0, 0),
                'human_clinical_abs_beta_ratio': (.928, .0005), 'human_clinical_SE_ratio': (1.019, .0005)}
    old_se = dict(zip(['whole_linear', 'whole_age_factor', 'whole_age_factor_slide',
                      'cortex_linear', 'cortex_age_factor', 'cortex_age_factor_slide'],
                     [1.572, 1.573, 1.800, 1.749, 1.648, 1.948]))
    for spec, value in precision.groupby('spec').SE_ratio.median().items():
        name = f'mouse_SE_ratio_{spec}'
        actual[name] = value
        expected[name] = (old_se[spec], .001)
    rows = [dict(metric=k, manuscript=v[0], computed=actual[k], tolerance=v[1],
                 matches=bool(np.isfinite(actual[k]) and abs(actual[k] - v[0]) <= v[1])) for k, v in expected.items()]
    pd.DataFrame(rows).to_csv(r / 'manuscript_comparison.csv', index=False)
    failures = dict(mouse_fits=int(status.failed.sum()), variance_fits=int(var.status.ne('available').sum()),
                    simulation_fits=int(sim.numerical_failed.sum()))
    unresolved = int(dream_diagnostics.gene_symbol.isna().sum())
    result = dict(status='completed_with_fit_failures' if any(failures.values()) else
                  ('completed_with_unresolved_warnings' if unresolved else 'completed'),
                  integrity_checks='passed', failures=failures,
                  dream_warning_records=len(dream_diagnostics), dream_unresolved_warnings=unresolved,
                  simulation_not_estimable=int(sim.not_estimable.sum()),
                  manuscript_anchor_differences=sum(not row['matches'] for row in rows),
                  historical_manuscript_reproduction_certified=False,
                  exact_SE_max_relative_error=float(error))
    dump(r / 'result_audit.json', result)
    return result
