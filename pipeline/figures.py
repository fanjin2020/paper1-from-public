"""Eight paper-aligned result figures, drawn only from computed CSV tables."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PALETTE = {'baseline': '#62696F', 'cage_fixed': '#C27435', 'cage_correlation': '#3B78A0',
           'dream_shared': '#58A99C', 'dream_native': '#258374',
           'mundlak_shared': '#BA85C1', 'mundlak_native': '#8B5793',
           'dream': '#258374', 'mundlak_dream': '#8B5793'}
LABELS = {'baseline': 'No cage', 'cage_fixed': 'Cage fixed', 'cage_correlation': 'Common corr.',
          'dream_shared': 'dream shared', 'dream_native': 'dream native',
          'mundlak_shared': 'Mundlak shared', 'mundlak_native': 'Mundlak native',
          'dream': 'dream', 'mundlak_dream': 'Mundlak'}


def clinical_pairs(frame):
    z = frame.loc[frame.contrast.eq('clinical_status')]
    return z.loc[z.model.eq('base')].merge(z.loc[z.model.eq('quality')],
        on=['gene_symbol', 'spec', 'context', 'contrast'], suffixes=('_base', '_quality'),
        validate='one_to_one')


def family_scope_labels(scopes):
    labels = {}
    for scope, g in scopes[scopes.model.eq('base') & scopes.denominator.eq('planned')].groupby('scope'):
        sizes = sorted(g.planned.dropna().astype(int).unique())
        size = f'{sizes[0]:,}' if len(sizes) == 1 else f'{sizes[0]:,}-{sizes[-1]:,}'
        labels[scope] = scope.replace('_', ' ')+'\n'+size+' / family'
    return labels


def highlighted_paired_design(frame):
    return (frame.n.eq(80) & frame.leverage.eq('ordinary') & frame.weights.eq('uniform')
            & np.isclose(frame.partial_correlation, .3))


def save(fig, frame, out, name):
    fig.savefig(out / (name + '.pdf'), bbox_inches='tight')
    fig.savefig(out / (name + '.png'), dpi=300, bbox_inches='tight')
    frame.to_csv(out / (name + '.csv'), index=False)
    plt.close(fig)


def errorbar(ax, x, z, metric, model, factor=1, mean=False):
    col = metric + '_mean' if mean else metric
    y = z[col].to_numpy() * factor
    lo, hi = metric + '_low', metric + '_high'
    err = None
    if lo in z and hi in z:
        err = np.maximum(0, np.vstack([z[col]-z[lo], z[hi]-z[col]]) * factor)
    elif metric == 'bias' and 'bias_MCSE' in z:
        err = 1.96*z.bias_MCSE.to_numpy()
    ax.errorbar(x, y, yerr=err, fmt='o', ms=4, capsize=2, color=PALETTE[model], label=LABELS[model])


def generate(core, matched, benchmark, out):
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
    c, m, b = core/'results', matched/'results', benchmark/'final_results'
    h = clinical_pairs(pd.read_csv(c/'human_results.csv.gz'))
    selected = h.loc[h.q_joint_base.lt(.05)]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), layout='constrained')
    axes[0].scatter(selected.q_joint_base, selected.q_joint_quality, s=20, color='#C27435')
    axes[0].axhline(.05, ls='--', color='grey')
    axes[0].set(xlabel='Base q', ylabel='Quality-adjusted q', title='A  Base-significant records')
    axes[1].scatter(h.beta_base, h.beta_quality, s=5, alpha=.3, color='#3B78A0')
    axes[1].scatter(selected.beta_base, selected.beta_quality, s=16, color='#C27435')
    vals = h[['beta_base','beta_quality']].to_numpy()
    lim = np.nanmax(np.abs(vals[np.isfinite(vals)]))
    axes[1].plot([-lim,lim],[-lim,lim], color='grey', ls='--')
    axes[1].set(xlabel='Base coefficient', ylabel='Adjusted coefficient', title='B  Clinical-status estimates')
    bins = [0,.01,.02,.03,.04,.05,.1,.2]
    counts = pd.DataFrame({arm: pd.cut(selected['q_joint_'+arm], bins, right=False, include_lowest=True).value_counts(sort=False)
                           for arm in ('base','quality')})
    counts.plot.bar(ax=axes[2], color=['#3B78A0','#C27435'], rot=0)
    axes[2].set_xticklabels(['<.01','.01-.02','.02-.03','.03-.04','.04-.05','.05-.10','.10-.20'], rotation=40, ha='right')
    axes[2].set(title='C  Selected-record q categories', xlabel='q', ylabel='Records')
    save(fig,h,out,'Figure1_human_estimates')

    family = pd.read_csv(c/'human_correction_scope_summary.csv')
    scope_labels = family_scope_labels(pd.read_csv(c/'correction_scope.csv'))
    fig, axes = plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for ax, metric in zip(axes, ['significant','original40_retained']):
        part = family.loc[family.denominator.eq('planned')] if metric=='significant' else family.loc[family.model.eq('quality')]
        column = 'model' if metric=='significant' else 'denominator'
        z = part.pivot(index='scope',columns=column,values=metric)
        z.plot.bar(ax=ax, rot=20, color=['#3B78A0','#C27435'])
        if metric == 'significant':
            ax.set_xticklabels([scope_labels[x] for x in z.index],rotation=15,ha='right')
        ax.set(xlabel='Testing-family definition',ylabel='Records',
               title='A  Clinical discoveries' if metric=='significant' else 'B  Original 40 retained')
    save(fig,family,out,'Figure2_testing_families')

    mouse = pd.read_csv(m/'mouse_method_summary.csv')
    comparison = pd.read_csv(m/'mouse_matched_method_summary.csv')
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    wide=mouse.pivot(index='spec',columns='model',values='significant_background')
    wide.plot.bar(ax=axes[0,0],rot=20,color=[PALETTE[x] for x in wide.columns])
    axes[0,0].set(ylabel='BH discoveries',xlabel='',title='A  Same animals and gene background')
    fits=pd.read_csv(m/'mouse_seven_arm_background_fits.csv.gz',usecols=['gene_symbol','spec','model','SE','status'])
    base=fits.loc[fits.model.eq('baseline')&fits.status.eq('available'),['gene_symbol','spec','SE']]
    se=fits.loc[fits.status.eq('available')&fits.model.ne('baseline')].merge(
        base,on=['gene_symbol','spec'],suffixes=('','_base'),validate='many_to_one')
    se['SE_ratio']=se.SE/se.SE_base
    ratios=se.groupby(['spec','model']).SE_ratio.median().reset_index()
    ratios.to_csv(out/'Figure3_SE_ratios.csv',index=False)
    for model,part in ratios.groupby('model'):
        axes[0,1].scatter([LABELS[model]]*len(part),part.SE_ratio,color=PALETTE[model],s=20)
    axes[0,1].tick_params(axis='x',rotation=30)
    axes[0,1].set(ylabel='Median paired SE ratio to no cage',title='B  Six specifications')
    axes[1,0].scatter(wide.cage_fixed,wide.mundlak_native,color=PALETTE['mundlak_native'])
    axes[1,0].set(xlabel='Fixed-effect discoveries',ylabel='Native Mundlak discoveries',title='C  Discovery counts')
    z=comparison.loc[comparison.model.eq('mundlak_native') & comparison.reference.isin(['cage_fixed','dream_native'])]
    for ref,part in z.groupby('reference'):
        axes[1,1].scatter([LABELS[ref]]*len(part),part.median_SE_ratio,label=LABELS[ref],color=PALETTE[ref])
    axes[1,1].axhline(1,color='grey',ls='--')
    axes[1,1].set(ylabel='Median paired Mundlak / comparator SE',title='D  Matched-gene contrasts')
    save(fig,mouse,out,'Figure3_mouse_models')
    comparison.to_csv(out/'Figure3_paired_comparisons.csv',index=False)

    from .publication_tables import native_precision, mouse_geometry
    precision=native_precision(pd.read_csv(c/'mouse_SE_decomposition_check.csv.gz'), pd.read_csv(c/'mouse_precision.csv.gz'))
    geometry=mouse_geometry(pd.read_csv(core/'prepared/mouse_units.csv'))
    columns=['SE_ratio','unscaled_ratio','posterior_SD_ratio']
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    axes=axes.ravel()
    for i,col in enumerate(columns):
        grouped=precision.groupby('spec')[col]
        mid=grouped.median();lo=grouped.quantile(.25);hi=grouped.quantile(.75)
        axes[0].errorbar(np.arange(len(mid))+(i-1.5)*.13,mid,yerr=[mid-lo,hi-mid],fmt='o',capsize=2,label=col)
    axes[0].set_xticks(np.arange(len(mid)),mid.index,rotation=25,ha='right')
    axes[0].axhline(1,color='grey',ls='--');axes[0].legend(fontsize=7)
    axes[0].set(ylabel='Median and gene IQR',title='A  Native-weight two-factor identity')
    groups=list(precision.groupby('spec',sort=True))
    axes[1].boxplot([g.reweight_factor for _,g in groups],tick_labels=[s for s,_ in groups],showfliers=False)
    axes[1].tick_params(axis='x',rotation=25)
    axes[1].axhline(1,color='grey',ls='--')
    axes[1].set(ylabel='Native/shared fixed-effect SE',title='B  Weight re-estimation')
    gg=geometry[geometry.spec.str.startswith('whole') & geometry.branch.isin(['baseline','cage_fixed'])]
    for ax,metric,title in zip(axes[2:],['genotype_residual_SS','residual_df'],['C  Residualized genotype information','D  Residual degrees of freedom']):
        gg.pivot(index='spec',columns='branch',values=metric).plot.bar(ax=ax,rot=20,color=[PALETTE['baseline'],PALETTE['cage_fixed']])
        ax.set(title=title,xlabel='',ylabel=metric.replace('_',' '))
    save(fig,pd.concat([precision.assign(record_type='gene'),geometry.assign(record_type='geometry')],ignore_index=True),out,'Figure4_precision_decomposition')

    paired=pd.read_csv(b/'paired_summary.csv')
    z=paired.loc[paired.status.eq('regular') & paired.mean_case.eq('quality_zero')].copy()
    z['highlighted_design']=highlighted_paired_design(z)
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for col,law in enumerate(['matched_gaussian','heteroskedastic_gaussian','matched_student_t5']):
        for row,metric in enumerate(['rejection','coverage']):
            ax=axes[row,col]
            for method,color in [('model_based','#3B78A0'),('HC3','#C27435')]:
                part=z.loc[z.error_law.eq(law)&z.method.eq(method)]
                ax.scatter(part.max_leverage,part[metric],s=17,label=method,color=color)
                marked=part[part.highlighted_design]
                ax.scatter(marked.max_leverage,marked[metric],s=85,marker='D',facecolors='none',
                           edgecolors='black',linewidths=1,
                           label='n=80, ordinary, equal; r=0.3' if method=='model_based' else '_nolegend_')
            ax.axhline(.05 if row==0 else .95,color='grey',ls='--')
            ax.set(xlabel='Maximum leverage',ylabel=metric,title=chr(65+row*3+col)+'  '+law.replace('_',' '))
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside upper center', ncol=3, fontsize=8, frameon=False)
    save(fig,z,out,'Figure5_paired_calibration')

    counts=pd.read_csv(m/'count_pipeline_summary.csv')
    z=counts.loc[counts.de_fraction.gt(0)]
    cases=z[['cage_log2_sd','cage_mean_amplitude']].drop_duplicates().sort_values(['cage_log2_sd','cage_mean_amplitude'])
    methods=['baseline','cage_correlation','dream_native','cage_fixed','mundlak_native']
    fig,axes=plt.subplots(3,len(cases),figsize=(14,8),layout='constrained',squeeze=False)
    for col,(sd,shift) in enumerate(cases.itertuples(index=False,name=None)):
        part=z.loc[z.cage_log2_sd.eq(sd)&z.cage_mean_amplitude.eq(shift)].set_index('model').reindex(methods)
        for row,metric in enumerate(['FDR','power','coverage']):
            ax=axes[row,col]
            for i,model in enumerate(methods):
                errorbar(ax,[i],part.loc[[model]],metric,model,mean=True)
            if metric!='power':ax.axhline(.05 if metric=='FDR' else .95,color='grey',ls='--')
            if metric=='FDR':
                value=part.loc['mundlak_native','FDR_mean']
                if np.isfinite(value) and value <= .10:
                    ax.annotate(f'{100*value:.2f}%',xy=(methods.index('mundlak_native'),value),
                                xytext=(-3,10),textcoords='offset points',ha='right',fontsize=8,
                                color=PALETTE['mundlak_native'])
                    ax.margins(y=.2)
            ax.set_xticks(range(len(methods)),[LABELS[x] for x in methods],rotation=50,ha='right')
            ax.set(ylabel=(chr(65+row)+'  '+metric) if col==0 else metric,
                   title=f'Cage SD {sd:g}; mean {shift:g}')
    save(fig,z,out,'Figure6_count_pipeline')

    cage=pd.read_csv(b/'cage_summary.csv')
    z=cage.loc[cage.layout.eq('actual_linear')&cage.truth.eq(0)]
    methods=['baseline','cage_correlation','cage_fixed','dream','mundlak_dream']
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for row,shift in enumerate([0,.5]):
        for col,metric in enumerate(['bias','rejection','mean_width']):
            ax=axes[row,col]
            for i,model in enumerate(methods):
                part=z.loc[z['shift'].eq(shift)&z.model.eq(model)].sort_values('rho_truth')
                errorbar(ax,part.rho_truth.to_numpy()+(i-2)*.009,part,metric,model)
            if col<2:ax.axhline(0 if col==0 else .05,color='grey',ls='--')
            ax.set(xlabel='Residual correlation',ylabel=metric,title=f'Cage mean {shift:g}')
            ax.set_xticks([0,.2])
    axes[0,0].legend(fontsize=7);save(fig,z,out,'Figure7_cage_benchmark')

    layouts=['n72_s6_m000_g36','n72_s6_m050_g36','n72_s6_m100_g36']
    z=cage.loc[cage.layout.isin(layouts)&cage.rho_truth.eq(.2)&cage['shift'].eq(0)]
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for ax,(metric,truth) in zip(axes.flat,[('rejection',0),('coverage',0),('rejection',.3),('mean_width',.3)]):
        for i,model in enumerate(methods):
            part=z.loc[z.truth.eq(truth)&z.model.eq(model)].set_index('layout').reindex(layouts)
            errorbar(ax,np.arange(3)+(i-2)*.07,part,metric,model)
        ax.set_xticks([0,1,2],['Pure','Half mixed','All mixed'])
        ax.set(ylabel=metric,title=f'True coefficient {truth:g}')
        if truth==0:ax.axhline(.05 if metric=='rejection' else .95,color='grey',ls='--')
    axes[0,0].legend(fontsize=7);save(fig,z,out,'Figure8_within_cage_allocation')
    pd.DataFrame([{'figure':p.stem,'pdf':p.name,'csv':p.with_suffix('.csv').name}
                  for p in sorted(out.glob('Figure*.pdf'))]).to_csv(out/'figure_index.csv',index=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('core','matched','benchmark','output'):
        p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    generate(args.core,args.matched,args.benchmark,args.output)


if __name__=='__main__':
    main()
