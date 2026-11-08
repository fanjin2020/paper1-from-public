"""Coverage and extreme-tail evidence from the completed paired benchmarks."""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA = Path(__file__).resolve().parent / 'data'
ARMS = ['estimated_shared', 'estimated_design', 'independent_shared',
        'independent_design', 'fixed_unit', 'gaussian_unit']
LABELS = ['Same counts\nshared weights', 'Same counts\nseparate weights',
          'Independent counts\nshared weights', 'Independent counts\nseparate weights',
          'Count response\nunit weights', 'Gaussian response\nunit weights']
CASES = ['s01', 's09', 's11', 's12', 's15', 's23']


def read(name):
    return pd.read_csv(DATA / ('paired_' + name + '.csv.gz'))


def arm_cloud(ax, z, metric, c, labels=True):
    for i, arm in enumerate(ARMS):
        for j, test in enumerate(['model_based', 'HC3']):
            group = z[z.arm.eq(arm) & z.test.eq(test)].set_index('scenario').loc[CASES]
            y = i + (j-.5)*.27 + np.linspace(-.075, .075, len(group))
            ax.scatter(group[metric]*100, y, color=c.PALETTE[test], marker=c.MARKERS[test],
                       s=24, edgecolors='white', linewidths=.4, zorder=3)
    ax.set_yticks(range(6), LABELS if labels else ['']*6)
    ax.set_ylim(5.6, -.6)
    ax.set_xlim(0, 16)
    ax.set_xticks([0, 5, 10, 15])
    ax.tick_params(axis='y', length=0, labelsize=8)
    c.ref(ax, 5, vertical=True)
    c.tidy(ax, 'x')


def fig10(c):
    z = read('tail_summary')
    fig = plt.figure(figsize=(c.WIDTH, 6.65))
    axes = [fig.add_axes(rect) for rect in [(.235,.53,.30,.37), (.665,.53,.30,.37),
                                          (.12,.12,.36,.27), (.615,.12,.35,.27)]]
    arm_cloud(axes[0], z, 'fdp_planned_050', c)
    arm_cloud(axes[1], z, 'bonferroni_fwer', c, labels=False)
    axes[0].set_xlabel('BH FDR (%)')
    axes[1].set_xlabel('Bonferroni FWER (%)')
    g = z[z.arm.eq('gaussian_unit')]
    for test in ['model_based', 'HC3']:
        a = g[g.test.eq(test)]
        axes[2].scatter(a.coverage95*100, a.fdp_planned_050*100, color=c.PALETTE[test],
                        marker=c.MARKERS[test], s=35, edgecolors='white', linewidths=.5)
    axes[2].set(xlabel='Gaussian 95% coverage (%)', ylabel='BH FDR (%)',
                xlim=(94.94,95.8), ylim=(0,14))
    axes[2].set_xticks([95,95.25,95.5,95.75], ['95.00','95.25','95.50','95.75'])
    axes[2].tick_params(axis='x', labelsize=7.8)
    c.ref(axes[2],95,vertical=True); c.ref(axes[2],5); c.tidy(axes[2])
    tail = read('tail_rates')
    tail = tail[tail.arm.eq('gaussian_unit') & tail.label.eq('bh_first')]
    for j,test in enumerate(['model_based','HC3']):
        a = tail[tail.test.eq(test)].set_index('scenario').loc[CASES]
        axes[3].scatter(np.arange(6)+(j-.5)*.22, a.ratio_to_uniform,
                        color=c.PALETTE[test], marker=c.MARKERS[test], s=32,
                        edgecolors='white', linewidths=.5, zorder=3)
    axes[3].set_xticks(range(6), CASES, rotation=45, ha='right')
    axes[3].set(xlabel='Gaussian control scenario', ylabel='Extreme-tail / uniform rate', ylim=(0,4))
    c.ref(axes[3],1); c.tidy(axes[3])
    for ax,letter in zip(axes,'abcd'): c.panel(ax,letter)
    c.legend(fig,['model_based','HC3'],y=.99)
    fig.text(.5,.025,'500 datasets per scenario  |  3,000 planned genes  |  zero paired differences',
             ha='center',fontsize=8)
    c.export(fig,'Fig10')


def figs14(c):
    z = read('count_summary')
    fig,axes=plt.subplots(2,2,figsize=(c.WIDTH,5.8))
    fig.subplots_adjust(left=.12,right=.98,bottom=.13,top=.88,wspace=.36,hspace=.44)
    for col,mode in enumerate(['shared','design_specific']):
        a=z[z.weight_mode.eq(mode)]
        for test in ['model_based','HC3']:
            b=a[a.test.eq(test)]
            axes[0,col].scatter(b.coverage95*100,b.fdp_planned_050*100,
                               color=c.PALETTE[test],marker=c.MARKERS[test],s=23,
                               alpha=.85,edgecolors='white',linewidths=.4)
            axes[1,col].scatter(b.fdp_planned_050*100,b.fdp_retained_050*100,
                               color=c.PALETTE[test],marker=c.MARKERS[test],s=23,
                               alpha=.85,edgecolors='white',linewidths=.4)
        axes[0,col].set(xlabel='95% interval coverage (%)',ylabel='Planned-family FDR (%)',
                        xlim=(94.2,96),ylim=(0,18))
        axes[0,col].set_title(['Shared weights','Separate weights'][col],pad=12)
        axes[1,col].set(xlabel='Planned-family FDR (%)',ylabel='Retained-family FDR (%)',
                        xlim=(0,18),ylim=(0,18))
        axes[1,col].plot([0,18],[0,18],color=c.GRAY,lw=.8,ls=':')
        for row in range(2):
            c.ref(axes[row,col],5)
            c.ref(axes[row,col],95 if row==0 else 5,vertical=True)
            c.tidy(axes[row,col])
            c.panel(axes[row,col],'abcd'[2*row+col])
    c.legend(fig,['model_based','HC3'],y=.99)
    fig.text(.5,.025,'24 scenarios  |  1,000 datasets per scenario  |  3,000 planned genes',
             ha='center',fontsize=8)
    c.export(fig,'FigS14')


def figs15(c):
    z=read('tail_summary')
    fig,axes=plt.subplots(2,3,figsize=(c.WIDTH,5.9),sharex=True,sharey=True)
    fig.subplots_adjust(left=.25,right=.98,bottom=.12,top=.88,wspace=.20,hspace=.37)
    for index,(ax,scenario) in enumerate(zip(axes.flat,CASES)):
        a=z[z.scenario.eq(scenario)]
        for i,arm in enumerate(ARMS):
            for j,test in enumerate(['model_based','HC3']):
                r=a[a.arm.eq(arm)&a.test.eq(test)].iloc[0]
                err=np.array([[r.fdp_planned_050-r.fdp_planned_050_ci_low],
                              [r.fdp_planned_050_ci_high-r.fdp_planned_050]])*100
                ax.errorbar([r.fdp_planned_050*100],[i+(j-.5)*.25],xerr=err,
                            color=c.PALETTE[test],marker=c.MARKERS[test],ms=3.6,
                            ls='none',capsize=1.6,elinewidth=.8)
        ax.set_yticks(range(6),LABELS,fontsize=7.8)
        ax.tick_params(axis='y',length=0)
        ax.set_ylim(5.55,-.55);ax.set_xlim(0,20);ax.set_xticks([0,5,10,15,20])
        ax.tick_params(axis='x',labelsize=7.5)
        r=a.iloc[0]
        ax.set_title(f'{scenario}: n = {r.n}',fontsize=9,pad=10)
        c.ref(ax,5,vertical=True);c.tidy(ax,'x');c.panel(ax,'abcdef'[index])
    c.legend(fig,['model_based','HC3'],y=.99)
    fig.text(.60,.035,'Planned-family BH FDR (%) with 95% Monte Carlo intervals',ha='center',fontsize=8)
    c.export(fig,'FigS15')
