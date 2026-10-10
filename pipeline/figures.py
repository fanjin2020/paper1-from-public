"""仅从本次统计输出作图。每幅图附带精简CSV，不依赖手填的论文数值。"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .io import require
from .audit import clinical_pairs

COLORS = ['#246B8E', '#B75939', '#39816E', '#8E6097', '#B08C2D']


def save(fig, frame, out, name):
    """同时输出矢量PDF、600dpi PNG和图示数据。"""
    fig.savefig(out / (name + '.pdf'), bbox_inches='tight')
    fig.savefig(out / (name + '.png'), dpi=600, bbox_inches='tight')
    frame.to_csv(out / (name + '.csv'), index=False)
    plt.close(fig)


def ecdf(ax, values, label, color):
    """绘制有限观测的经验累积分布；缺失值不作为零值进入分布。"""
    x = np.sort(np.asarray(values, float))
    x = x[np.isfinite(x)]
    if len(x):
        ax.step(x, np.arange(1, len(x)+1)/len(x), where='post', label=label, color=color)


def generate(work):
    """输出八组分析图及其CSV；图号与当前Word的对应关系需按图名核对。"""
    results = work / 'results'
    out = work / 'figures'
    out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    h = pd.read_csv(results / 'human_results.csv.gz')
    pair = clinical_pairs(h)
    pair['threshold_loss'] = (pair.q_joint_base < .05) & (pair.q_joint_quality >= .05)
    pair['coefficient_ratio'] = np.where(abs(pair.beta_base) > 1e-12, abs(pair.beta_quality/pair.beta_base), np.nan)
    pair['SE_ratio'] = pair.SE_quality/pair.SE_base
    sig = pair.loc[pair.q_joint_base < .05].copy()
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), layout='constrained')
    axes[0].scatter(pair.beta_base, pair.beta_quality, s=8, alpha=.35, color=COLORS[0])
    vals = pair[['beta_base', 'beta_quality']].to_numpy()
    finite = vals[np.isfinite(vals)]
    if len(finite):
        lim = max(abs(finite.min()), abs(finite.max()))
        axes[0].plot([-lim, lim], [-lim, lim], '--', color='.5', lw=1)
    axes[0].set(xlabel='Base coefficient', ylabel='Quality-adjusted coefficient', title='A  Clinical-status endpoints')
    axes[1].scatter(sig.q_joint_base, sig.q_joint_quality, s=15, color=COLORS[1])
    axes[1].axhline(.05, color='.4', ls='--', lw=1)
    axes[1].set(xlabel='Base adjusted P', ylabel='Quality-adjusted P', title='B  Base-significant endpoints')
    for col, color in zip(['coefficient_ratio', 'SE_ratio'], COLORS):
        ecdf(axes[2], sig[col], col.replace('_', ' '), color)
    axes[2].axvline(1, color='.5', ls='--'); axes[2].legend(fontsize=8)
    axes[2].set(xlabel='Paired ratio', ylabel='Cumulative proportion', title='C  Effect and precision')
    save(fig, pair, out, 'Figure1_human_effect_precision')

    delta = pd.read_csv(results / 'human_paired_changes.csv.gz')
    merged = sig.merge(delta, on=['gene_symbol', 'spec'], how='left', validate='one_to_one')
    fig, ax = plt.subplots(figsize=(8, max(3, min(14, len(merged)*.16))), layout='constrained')
    z = merged.sort_values('delta').reset_index(drop=True)
    for i, row in z.iterrows():
        if np.isfinite(row.delta) and np.isfinite(row.low):
            ax.plot([row.low, row.high], [i, i], color=COLORS[0], lw=.7)
            ax.scatter(row.delta, i, s=10, color=COLORS[1] if row.q_joint < .05 else COLORS[0])
    ax.axvline(0, color='.4', ls='--'); ax.set(xlabel='Quality minus base coefficient (95% CI)',
                                           ylabel='Base-significant endpoint rank')
    save(fig, merged, out, 'Figure2_paired_changes')

    precision = pd.read_csv(results / 'mouse_precision.csv.gz')
    specs = list(precision.spec.unique())
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), layout='constrained')
    for ax, spec in zip(axes.flat, specs):
        z = precision.loc[precision.spec == spec]
        for col, color in zip(['SE_ratio', 'design_factor', 'posterior_factor', 'reweight_factor'], COLORS):
            ecdf(ax, z[col], col.replace('_', ' '), color)
        ax.axvline(1, color='.6', ls='--', lw=.7)
        ax.set(title=spec.replace('_', ' '), xlabel='Factor', ylabel='Cumulative proportion')
        ax.legend(fontsize=7)
    save(fig, precision, out, 'Figure3_mouse_precision')

    variance = pd.read_csv(results / 'mouse_variance.csv.gz')
    valid_var = variance.loc[variance.status.eq('available')]
    summary = valid_var.groupby(['spec', 'weight_policy']).agg(available_genes=('gene_symbol', 'size'), median_ICC=('ICC', 'median'),
                zero_rate=('at_zero', 'mean'), singular_rate=('singular', 'mean'), converged_rate=('converged', 'mean')).reset_index()
    planned = variance.groupby(['spec', 'weight_policy']).size().rename('planned_genes').reset_index()
    summary = planned.merge(summary, on=['spec', 'weight_policy'], how='left')
    summary['available_genes'] = summary.available_genes.fillna(0).astype(int)
    summary['failed_genes'] = summary.planned_genes - summary.available_genes
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    for ax, policy in zip(axes, ['shared', 'native']):
        for i, (spec, z) in enumerate(valid_var.loc[valid_var.weight_policy.eq(policy)].groupby('spec')):
            ecdf(ax, z.ICC, spec, plt.cm.tab10(i))
        ax.set(xlabel='Cage ICC relative to average residual variance', ylabel='Cumulative proportion', title=policy.capitalize() + ' weights')
        ax.legend(fontsize=7)
    summary.to_csv(out / 'cage_variance_summary.csv', index=False)
    save(fig, variance[['gene_symbol', 'spec', 'weight_policy', 'status', 'ICC', 'at_zero', 'singular', 'converged']],
         out, 'Figure4_cage_variance')

    sim = pd.read_csv(results / 'simulation_summary.csv')
    z = sim.loc[(sim.truth == 0) & (sim.rho_truth == 0)].copy()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout='constrained')
    for i, (model, group) in enumerate(z.groupby('model')):
        axes[0].scatter(group.bias, group.rejection, label=model, color=COLORS[i], alpha=.7)
        axes[1].scatter(group.width, group.coverage, label=model, color=COLORS[i], alpha=.7)
        axes[2].scatter(group.failure_rate, group.rejection, label=model, color=COLORS[i], alpha=.7)
    axes[0].set(xlabel='Bias', ylabel='Null rejection probability')
    axes[1].set(xlabel='Mean interval width', ylabel='Conditional coverage')
    axes[2].set(xlabel='Failure rate', ylabel='Null rejection probability')
    axes[0].axhline(.05, color='.5', ls='--'); axes[1].axhline(.95, color='.5', ls='--')
    axes[0].legend(fontsize=8)
    save(fig, sim, out, 'Figure5_design_simulation')

    cal = pd.read_csv(results / 'calibration' / 'calibration_summary.csv')
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout='constrained')
    for i, (method, group) in enumerate(cal.groupby('method')):
        z = group.loc[group.mean_case != 'delta_two_reference_SD']
        axes[0].scatter(z.empirical_SD, z.RMS_SE, s=8, alpha=.4, label=method, color=COLORS[i])
        ecdf(axes[1], z.rejection, method, COLORS[i])
        ecdf(axes[2], group.loc[group.mean_case == 'delta_two_reference_SD', 'rejection'], method, COLORS[i])
    axes[0].set(xlabel='Empirical sampling SD', ylabel='RMS estimated SE')
    axes[1].axvline(.05, color='.5', ls='--');axes[1].set(xlabel='Null rejection probability', ylabel='Cumulative proportion')
    axes[2].set(xlabel='Power at two reference SD', ylabel='Cumulative proportion')
    axes[0].legend();save(fig, cal, out, 'Figure6_paired_test_calibration')

    mouse = pd.read_csv(results / 'mouse_results.csv.gz')
    rows = []
    for (spec, model), frame in mouse.groupby(['spec', 'model']):
        rows.append(dict(spec=spec, model=model, planned=len(frame), available=int(frame.p.notna().sum()),
                         significant_within_model=int(frame.q_within_model.lt(.05).sum()),
                         significant_base_fixed=int(frame.q_base_fixed_joint.lt(.05).sum())
                         if model in ['baseline', 'cage_fixed'] else np.nan))
    counts = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout='constrained')
    wide = counts.pivot(index='spec', columns='model', values='significant_within_model')
    wide.plot.bar(ax=axes[0], rot=25, color=COLORS, width=.8)
    axes[0].set(ylabel='Genes passing within-model BH', xlabel='')
    axes[0].legend(fontsize=7)
    core = counts.loc[counts.model.isin(['baseline', 'cage_fixed'])]
    core.pivot(index='spec', columns='model', values='significant_base_fixed').plot.bar(
        ax=axes[1], rot=25, color=COLORS[:2])
    axes[1].set(ylabel='Genes passing joint two-model BH', xlabel='')
    for ax in axes:
        ax.tick_params(axis='x', labelsize=7)
    save(fig, counts, out, 'Figure7_mouse_testing_families')

    implementation = pd.read_csv(results / 'mouse_SE_implementation_audit.csv.gz')
    implementation['contrast_to_direct'] = implementation.SE_contrast_ordinary / implementation.SE_direct_ordinary
    implementation['robust_to_ordinary'] = implementation.SE_direct_robust / implementation.SE_direct_ordinary
    diagnostic = implementation.groupby(['spec', 'branch']).agg(
        genes=('gene_symbol', 'size'), median_contrast_ratio=('contrast_to_direct', 'median'),
        median_robust_ratio=('robust_to_ordinary', 'median')).reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    for i, (branch, frame) in enumerate(implementation.groupby('branch')):
        ecdf(axes[0], frame.contrast_to_direct, branch, COLORS[i])
        ecdf(axes[1], frame.robust_to_ordinary, branch, COLORS[i])
    for ax in axes:
        ax.axvline(1, color='.5', ls='--', lw=.8)
        ax.legend(fontsize=8)
        ax.set(ylabel='Cumulative proportion')
    axes[0].set(xlabel='Contrast-transform / direct-coefficient SE')
    axes[1].set(xlabel='Robust / ordinary moderation SE')
    diagnostic.to_csv(out / 'SE_implementation_summary.csv', index=False)
    save(fig, implementation[['gene_symbol', 'spec', 'branch', 'contrast_to_direct', 'robust_to_ordinary']],
         out, 'Figure8_SE_implementation')

    # 计数核对表报告实际得到的值，不把论文旧数字当成计算目标。
    rows = [dict(metric='frontal_panel_genes', value=len(pd.read_csv(results/'panel.csv'))),
            dict(metric='human_clinical_base_significant', value=int((pair.q_joint_base<.05).sum())),
            dict(metric='human_clinical_quality_significant', value=int((pair.q_joint_quality<.05).sum())),
            dict(metric='human_threshold_losses', value=int(pair.threshold_loss.sum())),
            dict(metric='human_selected_abs_beta_ratio_median', value=sig.coefficient_ratio.median()),
            dict(metric='human_selected_SE_ratio_median', value=sig.SE_ratio.median())]
    pd.DataFrame(rows).to_csv(out/'manuscript_numbers.csv', index=False)
    h.groupby(['contrast', 'model']).agg(planned=('p', 'size'), available=('p', 'count'),
        significant=('q_joint', lambda x: int(x.lt(.05).sum()))).reset_index().to_csv(
            out/'human_contrast_counts.csv', index=False)
    pd.DataFrame([dict(figure=p.stem, pdf=p.name, source_csv=p.with_suffix('.csv').name)
                  for p in sorted(out.glob('Figure*.pdf'))]).to_csv(out/'figure_index.csv', index=False)
    print('新结果图表已生成：', out, flush=True)
