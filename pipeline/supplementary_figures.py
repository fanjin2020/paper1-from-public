"""Twelve supplementary figures with stable manuscript numbering and sources."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .figures import PALETTE, LABELS, save, errorbar
from .publication_tables import ICC_LABELS


def forest(ax, z, labels, beta, low, high):
    good=np.isfinite(z[[beta,low,high]].to_numpy(float)).all(axis=1)
    a=z.loc[good]; y=np.flatnonzero(good)
    ax.hlines(y,a[low],a[high],color='#3B78A0',lw=1)
    ax.scatter(a[beta],y,s=15,color='#3B78A0')
    ax.set_yticks(range(len(z)),labels,fontsize=7)
    ax.axvline(0,color='grey',ls='--',lw=.7)
    ax.set_xlabel('Coefficient and pointwise 95% interval')


def generate(report):
    report=Path(report); out=report/'figures'; out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,
                         'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    read=lambda name: pd.read_csv(report/name,low_memory=False)
    index=[]
    def finish(fig,z,name,source,definition):
        save(fig,z,out,name)
        index.append(dict(figure=name,numeric_source=source,rows=len(z),definition=definition))

    src='S05/donors.csv';z=read(src)
    fig,axes=plt.subplots(1,3,figsize=(11,4),layout='constrained')
    for ax,field in zip(axes,['PMI','RIN','pH']):
        data=[z.loc[z.clinical_group.eq(k),field].dropna().to_numpy() for k in ['No_dementia','Dementia']]
        ax.boxplot(data,tick_labels=['No dementia','Dementia'],showfliers=False)
        for i,y in enumerate(data,1):ax.scatter(i+np.linspace(-.12,.12,len(y)),y,s=12,color=['#3B78A0','#C27435'][i-1],alpha=.7)
        ax.set_ylabel(field)
    finish(fig,z,'FigureS1_donor_quality',src,'One dot per donor; boxes 1.5-IQR whiskers.')

    src='S08/SEA_astrocyte_candidates.csv';z=read(src)
    fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
    forest(ax,z,z.gene_symbol+' / '+z.contrast+' / '+z.model,'beta','low','high')
    ax.collections[-1].remove()
    for model,marker in [('base','o'),('quality','s')]:
        good=z.model.eq(model)&z.beta.notna();ax.scatter(z.loc[good,'beta'],np.flatnonzero(good),marker=marker,s=18,label=model)
    ax.legend()
    finish(fig,z,'FigureS2_astrocyte_estimates',src,'Base and quality-adjusted fits, identical eligible donors.')

    src='S13/publication_mouse_consistency.csv';z=read(src)
    panel=read('S12/publication_three_mouse_arms.csv')
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    counts=z.significant.value_counts().sort_index();axes[0].bar(counts.index,counts,color='#3B78A0')
    axes[0].set(xlabel='Significant specifications out of 18',ylabel='Genes',title='A  Complete panel positions')
    focal=panel[panel.gene_symbol.isin(['GFAP','GSDMD','IGFBP5','EMP3'])].copy()
    focal['discovery']=focal.q_three_model_joint.lt(.05)
    a=focal.groupby('gene_symbol').discovery.sum();axes[1].bar(a.index,a,color='#C27435')
    axes[1].set(ylabel='Significant specifications out of 18',title='B  Focal genes')
    finish(fig,pd.concat([z.assign(record_type='all_genes'),focal.assign(record_type='focal')],ignore_index=True),'FigureS3_mouse_consistency',src,'Three-arm joint planned family; focal values also in S12.')

    src='S15/ICC_precision.csv';z=read(src);strata=read('S15/ICC_strata.csv')
    specs=sorted(z.spec.unique());fig,axes=plt.subplots(2,len(specs),figsize=(18,7),layout='constrained',squeeze=False)
    for col,spec in enumerate(specs):
        a=z[z.spec.eq(spec)&z.weight_policy.eq('shared')].dropna(subset=['ICC','SE_ratio'])
        axes[0,col].hexbin(a.ICC,a.SE_ratio,gridsize=35,mincnt=1,cmap='Blues')
        axes[0,col].set(title=spec,xlabel='Conditional ICC',ylabel='Shared dream/base SE')
        for policy,color in [('shared','#3B78A0'),('native','#258374')]:
            g=strata[strata.spec.eq(spec)&strata.weight_policy.eq(policy)].set_index('ICC_bin').reindex(ICC_LABELS)
            axes[1,col].plot(range(4),g.median_SE_ratio,'o-',label=policy,color=color)
        axes[1,col].set_xticks(range(4),ICC_LABELS,rotation=35,ha='right')
        axes[1,col].axhline(1,ls='--',color='grey');axes[1,col].set_ylabel('Median dream/base SE')
    axes[1,0].legend(fontsize=7)
    finish(fig,z,'FigureS4_ICC_precision',src,'dream/base by matching weight policy; ICC cut points .01, .10, .20.')

    src='S17/cage_summary.csv';z=read(src);z=z[z.rho_truth.eq(.2)]
    layouts=['actual_linear','n72_s6_m000_g36','n72_s6_m050_g36','n72_s6_m100_g36']
    methods=['baseline','cage_correlation','cage_fixed','dream','mundlak_dream']
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    for row,truth in enumerate([0,.3]):
        for col,shift in enumerate([0,.5]):
            ax=axes[row,col]
            for i,model in enumerate(methods):
                a=z[z.truth.eq(truth)&z['shift'].eq(shift)&z.model.eq(model)].set_index('layout').reindex(layouts)
                errorbar(ax,np.arange(4)+(i-2)*.09,a,'rejection',model)
            ax.set_xticks(range(4),['Observed','Pure','Half mixed','All mixed'])
            ax.set(ylabel='Rejection probability',title=f'True effect {truth:g}; mean shift {shift:g}')
            if truth==0:ax.axhline(.05,color='grey',ls='--')
    axes[0,0].legend(fontsize=7)
    finish(fig,z,'FigureS5_cage_rejection',src,'Five-model comparison at residual correlation .20, with Wilson intervals.')

    src='S04/publication_continuous.csv';z=read(src)
    z=z[z.collection.eq('human_17064')&z.model.eq('clinical_status')].copy()
    contexts=sorted(z.context.unique());fig,axes=plt.subplots(1,3,figsize=(14,8),layout='constrained')
    for ax,field,title in zip(axes[:2],['log_abs_beta_ratio','minus_log_SE_ratio'],['Log coefficient change','Negative log SE change']):
        a=z.groupby('context')[field];mid=a.median().reindex(contexts);lo=a.quantile(.25).reindex(contexts);hi=a.quantile(.75).reindex(contexts)
        ax.errorbar(mid,range(len(mid)),xerr=[mid-lo,hi-mid],fmt='o',capsize=2,color='#3B78A0')
        ax.set_yticks(range(len(mid)),contexts,fontsize=7);ax.axvline(0,color='grey',ls='--');ax.set_title(title)
    a=z.groupby('context').agg(available=('beta_base','count'),selected=('original_q_base',lambda q:q.lt(.05).sum())).reindex(contexts)
    a.plot.barh(ax=axes[2],color=['#3B78A0','#C27435']);axes[2].set_title('Available and selected')
    finish(fig,z,'FigureS6_human_continuous_changes',src,'All decomposable clinical records, not only selected discoveries.')

    src='S16/publication_gene_design_errors.csv';z=read(src);d=read('S16/publication_design_display.csv')
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for r in d.itertuples():
        axes[0].errorbar(r.predicted,r.observed_median,yerr=[[r.observed_median-r.observed_q25],[r.observed_q75-r.observed_median]],fmt='o',label=r.spec)
    lim=[min(d.predicted.min(),d.observed_q25.min())*.95,max(d.predicted.max(),d.observed_q75.max())*1.05]
    axes[0].plot(lim,lim,'--',color='grey');axes[0].set(xlabel='Design-only prediction',ylabel='Observed SE median and IQR');axes[0].legend(fontsize=6)
    groups=list(z.groupby('spec'));axes[1].boxplot([a.design_only_log_error for _,a in groups],tick_labels=[s for s,_ in groups],whis=(5,95),showfliers=False)
    axes[1].tick_params(axis='x',rotation=30);axes[1].axhline(0,ls='--',color='grey');axes[1].set_ylabel('Log observed/predicted SE ratio')
    finish(fig,z,'FigureS7_design_precision',src,'Expression-free prediction against background genes; whiskers 5th-95th percentiles.')

    src='S20/donor_conditioned_calibration.csv';cal=read(src)
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    z=cal[cal.error_law.eq('matched_gaussian') & cal.mean_case.eq('quality_nonzero_delta_zero')]
    for i,method in enumerate(['model_based','HC3']):
        a=z[z.method.eq(method)];axes[0].scatter(np.full(len(a),i+1),a.rejection,s=8,alpha=.4,label=method)
    a=z[z.method.eq('model_based')];axes[0].scatter(np.zeros(len(a)),a.joint_F_rejection,s=8,alpha=.4)
    axes[0].set_xticks([0,1,2],['Joint F','Paired model','Paired HC3']);axes[0].axhline(.05,ls='--',color='grey');axes[0].set_ylabel('Rejection probability')
    a=cal[cal.mean_case.eq('quality_zero')&cal.method.eq('HC3')]
    for law,g in a.groupby('error_law'):axes[1].scatter(g.max_leverage,g.RMS_SE/g.empirical_SD,s=8,alpha=.5,label=law)
    axes[1].axhline(1,ls='--',color='grey');axes[1].set(xlabel='Maximum leverage',ylabel='HC3 RMS SE / empirical SD');axes[1].legend(fontsize=6)
    finish(fig,cal,'FigureS8_joint_and_paired_tests',src,'Cancellation mean for joint F; zero quality effects for HC3 leverage diagnostic.')

    src='S04/publication_selected_log_components.csv';selected=read(src)
    fig,ax=plt.subplots(figsize=(11,6),layout='constrained');ax.axis('off')
    steps=[('Match the comparison','Biological units, target coding, response preprocessing'),
           ('Check coefficient change','Paired contrast with shared-response covariance'),
           ('Separate precision components','Weighted comparison, residual variance, weight re-estimation'),
           ('Check grouped-design information','Within-cage overlap, rank and residual degrees of freedom'),
           ('Report decisions and availability','Declared BH family; tested, unavailable and selected positions')]
    for i,(title,body) in enumerate(steps):
        y=.9-i*.19
        ax.text(.5,y,title+'\n'+body,ha='center',va='center',bbox=dict(boxstyle='square,pad=.5',facecolor='#F0F5F8',edgecolor='#3B78A0'),fontsize=10)
        if i<4:ax.annotate('',xy=(.5,y-.13),xytext=(.5,y-.06),arrowprops=dict(arrowstyle='->',color='#3B78A0'))
    finish(fig,pd.DataFrame(steps,columns=['step','purpose']),'FigureS9_diagnostic_workflow',src,'Protocol diagram, not an additional estimate; illustrative input source listed.')

    z=selected.sort_values('log_abs_t_ratio');fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    axes[0].scatter(z.log_abs_beta_ratio,z.minus_log_SE_ratio,s=20,color='#3B78A0')
    lim=np.nanmax(np.abs(z[['log_abs_beta_ratio','minus_log_SE_ratio']].to_numpy()))*1.1
    axes[0].plot([-lim,lim],[-lim,lim],'--',color='grey');axes[0].plot([-lim,lim],[lim,-lim],'--',color='grey')
    axes[0].set(xlabel='Log absolute coefficient ratio',ylabel='Negative log SE ratio')
    for field,label in [('log_abs_beta_ratio','Coefficient'),('minus_log_SE_ratio','Precision'),('log_abs_t_ratio','Sum')]:axes[1].plot(range(len(z)),z[field],'.-',ms=3,label=label)
    axes[1].set(xlabel='Selected record, ordered by log t change',ylabel='Log contribution');axes[1].legend()
    finish(fig,z,'FigureS10_selected_log_components',src,'Descriptive decomposition; not a coefficient-equality test.')

    src='S20/donor_conditioned_calibration.csv';z=cal.copy();z['SE_ratio']=z.RMS_SE/z.empirical_SD
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for ax,(metric,case,title,ref) in zip(axes.flat,[('rejection','quality_zero','A  Null rejection',.05),('coverage','quality_zero','B  Coverage',.95),('SE_ratio','quality_zero','C  SE calibration',1),('rejection','delta_two_reference_SD','D  Two-SD power',None)]):
        g=list(z[z.mean_case.eq(case)&z.method.isin(['model_based','HC3'])].groupby(['error_law','method']))
        ax.boxplot([a[metric] for _,a in g],tick_labels=['\n'.join(k) for k,_ in g],whis=(0,100),showfliers=False)
        for i,(_,a) in enumerate(g,1):ax.scatter(i+np.linspace(-.12,.12,len(a)),a[metric],s=3,alpha=.25,color='#3B78A0')
        ax.tick_params(axis='x',labelsize=6,rotation=20);ax.set_title(title)
        if ref is not None:ax.axhline(ref,ls='--',color='grey')
    finish(fig,z,'FigureS11_donor_calibration',src,'Donor-conditioned designs; A-C null, D two-reference-SD alternative.')

    src='S20/general_calibration.csv';z=read(src);z=z[z.status.eq('regular')&z.error_law.eq('matched_gaussian')]
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    cancel=z[z.mean_case.eq('quality_nonzero_delta_zero')]
    for _,g in cancel.groupby('design'):
        a=g.set_index('method');axes[0].plot([0,1,2],[a.loc['model_based','joint_F_rejection'],a.loc['model_based','rejection'],a.loc['HC3','rejection']],'.-',alpha=.3,color='#3B78A0')
    axes[0].set_xticks([0,1,2],['Joint F','Paired model','Paired HC3']);axes[0].axhline(.05,ls='--',color='grey');axes[0].set_ylabel('Rejection probability')
    for i,case in enumerate(['delta_one_reference_SD','delta_two_reference_SD']):
        for j,method in enumerate(['model_based','HC3']):
            a=z[z.mean_case.eq(case)&z.method.eq(method)].rejection;x=i*3+j
            axes[1].scatter(x+np.linspace(-.12,.12,len(a)),a,s=12,color=['#3B78A0','#C27435'][j]);axes[1].hlines(a.median(),x-.2,x+.2,color='black')
    axes[1].set_xticks([0,1,3,4],['1 SD\nModel','1 SD\nHC3','2 SD\nModel','2 SD\nHC3']);axes[1].set_ylabel('Pointwise power')
    finish(fig,z,'FigureS12_cancellation_power',src,'General fixed-design grid, separate from donor-conditioned simulations.')
    pd.DataFrame(index).to_csv(out/'figure_index.csv',index=False)
    return index
