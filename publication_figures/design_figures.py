"""Publication figures for the completed design-geometry benchmark."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


def fig9(p):
    designs = p.read('extension_design_diagnostics')
    counts = p.read('extension_counts_summary')
    hold = p.read('extension_heldout_rule_evaluation')
    empirical = p.read('extension_empirical_design_diagnostics').iloc[0]
    fig, axs = plt.subplots(2, 2, figsize=(p.WIDTH, 7.6))
    fig.subplots_adjust(left=.11, right=.98, bottom=.14, top=.88,
                        wspace=.40, hspace=.52)
    ax = axs[0, 0]
    for cages, color, marker in [(12,p.BLUE,'o'),(18,p.TEAL,'s'),(36,p.PURPLE,'^')]:
        z=designs[designs.groups.eq(cages)].sort_values('ewi')
        ax.plot(z.ewi,z.conditional_ewi,marker=marker,color=color,ms=5,
                lw=1,label=f'{cages} cages')
    ax.scatter([empirical.ewi],[empirical.conditional_ewi],marker='D',
               color=p.ORANGE,s=35,zorder=4,label='Observed layout')
    ax.plot([0,1],[0,1],ls=':',color=p.GRAY,lw=.8,zorder=0)
    ax.set(xlabel='Unadjusted within-group information',
           ylabel='Covariate-adjusted information ratio',xlim=(-.04,1.04),ylim=(-.04,1.04))
    ax.legend(loc='upper left',fontsize=8.2,handlelength=1.2)
    p.tidy(ax);p.panel(ax,'a')

    ax=axs[0,1]
    for i,m in enumerate(p.COUNT_MODELS):
        z=counts[counts.design.eq('k18_w1') & counts.model.eq(m)].set_index('scenario')
        values=100*z.loc[[4,5],'empirical_fdr'].to_numpy(float)
        ax.plot([0,1],values,color=p.PALETTE[m],marker=p.MARKERS[m],
                ms=5,lw=1.2)
    p.ref(ax,5)
    ax.set_xticks([0,1],['No cage-mean shift','Cage-mean shift'])
    ax.tick_params(axis='x',labelsize=8.2)
    ax.set(ylabel='False discovery rate (%)',ylim=(0,72),xlim=(-.17,1.17))
    ax.text(.03,.95,'18 cages; information ratio 0.333\n10% nonnull; latent cage SD 0',
            transform=ax.transAxes,va='top',fontsize=8.2)
    p.tidy(ax);p.panel(ax,'b')

    ax=axs[1,0]
    designs_order=[f'k18_w{i}' for i in range(4)]
    status={'within_effect_not_estimable':0,'no_supported_selection':1,'conditionally_eligible':2}
    grid=np.array([[status[hold[hold.design.eq(d)&hold.scenario.eq(s)].iloc[0].status]
                    for s in range(8)] for d in designs_order])
    cmap=ListedColormap(['#C5CDD3','#F2D7B9','#9ACBC3'])
    ax.imshow(grid,cmap=cmap,norm=BoundaryNorm([-.5,.5,1.5,2.5],3),aspect='auto')
    short={'baseline':'B','cage_correlation':'C','cage_fixed':'F','dream_native':'D','mundlak_native':'M'}
    for i,d in enumerate(designs_order):
        for s in range(8):
            r=hold[hold.design.eq(d)&hold.scenario.eq(s)].iloc[0]
            label=short[r.model] if r.status=='conditionally_eligible' else ('NE' if grid[i,s]==0 else 'A')
            ax.text(s,i,label,ha='center',va='center',fontsize=8.2,color=p.INK)
    ax.set_xticks(range(8),[str(s+1) for s in range(8)])
    ax.set_yticks(range(4),['0','1/3','2/3','1'])
    ax.set(xlabel='Generating scenario',ylabel='Unadjusted information ratio')
    ax.set_xticks(np.arange(-.5,8,1),minor=True)
    ax.set_yticks(np.arange(-.5,4,1),minor=True)
    ax.grid(which='minor',color='white',lw=1.2)
    ax.tick_params(which='minor',length=0)
    fig.legend(handles=[Patch(facecolor=cmap(i),label=l) for i,l in enumerate(
        ['NE: not estimable','A: abstention','Eligible'])],loc='upper center',
        bbox_to_anchor=(.50,.045),ncol=3,fontsize=8.2)
    p.panel(ax,'c')

    ax=axs[1,1]
    selected=hold[hold.status.eq('conditionally_eligible')].sort_values(['design','scenario'])
    labels=[]
    for i,r in enumerate(selected.itertuples()):
        ax.errorbar(100*r.heldout_fdr,i,xerr=np.array([[r.heldout_fdr-r.heldout_fdr_ci_lower],
                    [r.heldout_fdr_ci_upper-r.heldout_fdr]])*100,
                    fmt=p.MARKERS[r.model],color=p.PALETTE[r.model],ms=4,
                    elinewidth=.9,capsize=1.7)
        labels.append(f'{r.design[-1]}/3 | {r.scenario+1}')
    ax.set_yticks(range(len(selected)),labels,fontsize=8.0)
    ax.invert_yaxis()
    ax.set(xlabel='Held-out FDR (%)',ylabel='Information level | scenario',xlim=(-.5,18))
    p.ref(ax,5,vertical=True)
    ax.axvline(7.5,color=p.ORANGE,ls=':',lw=1)
    p.tidy(ax,'x');p.panel(ax,'d')
    handles=[Line2D([],[],color=p.PALETTE[m],marker=p.MARKERS[m],lw=1,
             label=p.LABELS[m].replace(' native','')) for m in p.COUNT_MODELS]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.53,1.0),
               ncol=3,fontsize=8)
    p.export(fig,'Fig9')


def figs13(p):
    z=p.read('extension_gaussian_summary')
    designs=p.read('extension_design_diagnostics')
    z=z.merge(designs[['design','ewi']],on='design',validate='many_to_one')
    models=['baseline','re_known_covariance','fixed','mundlak_known_covariance']
    labels=['No cage','Known-covariance RE','Cage fixed','Known-covariance Mundlak']
    colors=[p.GRAY,p.TEAL,p.ORANGE,p.PURPLE]
    markers=['o','^','D','P']
    fig,axs=plt.subplots(1,2,figsize=(p.WIDTH,3.65))
    fig.subplots_adjust(left=.11,right=.98,bottom=.20,top=.77,wspace=.40)
    for m,color,marker,label in zip(models,colors,markers,labels):
        a=z[z.model.eq(m)].dropna(subset=['empirical_positive_gamma_bias','analytic_positive_gamma_bias'])
        axs[0].scatter(a.analytic_positive_gamma_bias,a.empirical_positive_gamma_bias,
                      marker=marker,color=color,s=20,alpha=.7,linewidths=.3,edgecolors='white',label=label)
        a=a[a.scenario.eq(5)].sort_values(['ewi','design'])
        axs[1].scatter(a.ewi,a.empirical_positive_gamma_bias,marker=marker,
                      color=color,s=25,alpha=.8,linewidths=.3,edgecolors='white')
    axs[0].plot([-.02,.52],[-.02,.52],color=p.INK,ls='--',lw=.8)
    axs[0].set(xlabel='Analytical coefficient bias',ylabel='Mean simulated coefficient bias',
               xlim=(-.03,.53),ylim=(-.03,.53))
    axs[1].set(xlabel='Unadjusted information ratio',ylabel='Mean simulated coefficient bias',
               xlim=(-.04,1.04),ylim=(-.03,.53))
    p.ref(axs[1],0)
    for i,ax in enumerate(axs):
        p.panel(ax,'ab'[i]);p.tidy(ax)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',
               bbox_to_anchor=(.52,.99),ncol=2,fontsize=8)
    p.export(fig,'FigS13')
