"""Reproduce the manuscript's ten main and fifteen supplementary figures.

All numerical inputs are fixed source-table extracts. This script performs only
display aggregation; it does not refit models or rerun simulations.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.ticker import MaxNLocator
from matplotlib.text import Text

ROOT = Path(__file__).resolve().parent
WIDTH = 174 / 25.4
BLUE, ORANGE, TEAL, PURPLE = '#3976AF', '#D48642', '#198F82', '#885BA3'
INK, GRAY, GRID = '#263442', '#758391', '#E7EDF1'
PALETTE = {
    'baseline': GRAY, 'cage_correlation': BLUE, 'cage_fixed': ORANGE,
    'dream': TEAL, 'dream_native': TEAL, 'dream_shared': '#79BCAD',
    'mundlak_dream': PURPLE, 'mundlak_native': PURPLE, 'mundlak_shared': '#BE9ACB',
    'model_based': BLUE, 'HC3': ORANGE,
}
MARKERS = {'baseline':'o','cage_correlation':'s','cage_fixed':'D',
           'dream':'^','dream_native':'^','dream_shared':'v',
           'mundlak_dream':'P','mundlak_native':'P','mundlak_shared':'X',
           'model_based':'o','HC3':'s'}
LABELS = {'baseline':'No cage','cage_correlation':'Common correlation','cage_fixed':'Cage fixed',
          'dream':'dream','dream_native':'dream native','dream_shared':'dream shared',
          'mundlak_dream':'Mundlak','mundlak_native':'Mundlak native','mundlak_shared':'Mundlak shared',
          'model_based':'Model-based','HC3':'HC3'}
MODELS = ['baseline','cage_correlation','dream','cage_fixed','mundlak_dream']
COUNT_MODELS = ['baseline','cage_correlation','dream_native','cage_fixed','mundlak_native']
ARMS = ['baseline','cage_correlation','dream_shared','dream_native','cage_fixed','mundlak_shared','mundlak_native']
SPECS = ['whole_linear','whole_age_factor','whole_age_factor_slide',
         'cortex_linear','cortex_age_factor','cortex_age_factor_slide']
SPEC_LABELS = ['Whole\nLinear age','Whole\nAge factor','Whole\nAge + slide',
               'Cortex\nLinear age','Cortex\nAge factor','Cortex\nAge + slide']
SPEC_SHORT = dict(zip(SPECS,['Whole: linear','Whole: age factor','Whole: age + slide',
                           'Cortex: linear','Cortex: age factor','Cortex: age + slide']))
SPEC_COLORS = [BLUE,TEAL,'#81A794',PURPLE,ORANGE,'#AE6D7D']
LAWS = ['matched_gaussian','heteroskedastic_gaussian','matched_student_t5']
LAW_LABELS = ['Matched\nGaussian','Heteroskedastic\nGaussian','Matched\nStudent t5']
LAYOUTS = ['n72_s6_m000_g36','n72_s6_m050_g36','n72_s6_m100_g36']
OUT = None
TIFF = False
QA = []


def style():
    plt.rcParams.update({
        'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':8.7,'axes.labelsize':8.7,'xtick.labelsize':8.2,'ytick.labelsize':8.2,
        'axes.titlesize':9.5,'legend.fontsize':8.3,'text.color':INK,'axes.labelcolor':INK,
        'xtick.color':INK,'ytick.color':INK,'axes.edgecolor':'#97A3AD',
        'axes.linewidth':.65,'axes.spines.top':False,'axes.spines.right':False,
        'axes.axisbelow':True,'xtick.major.size':3,'ytick.major.size':3,
        'grid.color':GRID,'grid.linewidth':.55,'lines.linewidth':1.2,
        'legend.frameon':False,'legend.handlelength':1.3,'legend.columnspacing':1.25,
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
        'savefig.facecolor':'white','figure.facecolor':'white','axes.facecolor':'white',
    })


def read(name):
    return pd.read_csv(ROOT/'data'/(name+'.csv.gz'),low_memory=False)


def verify_inputs():
    manifest=pd.read_csv(ROOT/'source_index.csv')
    for row in manifest.itertuples():
        p=ROOT/row.file
        if hashlib.sha256(p.read_bytes()).hexdigest()!=row.sha256:
            raise ValueError('Input checksum mismatch: '+row.file)
    h=read('clinical_pairs').dropna()
    selected=h[h.q_joint_base.lt(.05)]
    bins=[0,.01,.02,.03,.04,.05,.10,.20]
    assert len(h)==3744 and len(selected)==40
    assert np.histogram(selected.q_joint_base,bins)[0].tolist()==[0,2,7,15,16,0,0]
    assert np.histogram(selected.q_joint_quality,bins)[0].tolist()==[0,0,0,0,0,15,25]
    c=read('count_summary')
    z=c[c.cage_log2_sd.eq(0)&c.cage_mean_amplitude.eq(.5)&c.de_fraction.gt(0)].set_index('model')
    for model, expected in [('dream_native',57.27),('cage_fixed',4.54),('mundlak_native',5.69)]:
        assert round(z.loc[model,'FDR_mean']*100,2)==expected


def panel(ax, letter, title=None):
    offset = .035 / ax.get_position().width
    ax.text(-offset,1.06,letter,transform=ax.transAxes,fontweight='bold',fontsize=11,
            va='bottom',ha='left',clip_on=False)
    if title:
        ax.set_title(title,loc='left',pad=12,fontweight='semibold')


def tidy(ax, axis='y'):
    ax.grid(axis=axis)
    ax.tick_params(pad=3)
    ax.margins(x=.08)


def ref(ax, value, vertical=False):
    f=ax.axvline if vertical else ax.axhline
    f(value,color='#7C8791',ls=(0,(3,2)),lw=.85,zorder=1)


def legend(fig, methods, y=.98, ncol=None):
    handles=[Line2D([],[],color=PALETTE[m],marker=MARKERS[m],ls='',ms=5,label=LABELS[m]) for m in methods]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,y),
               ncol=ncol or len(methods))


def points(ax,x,z,metric,model,scale=1,mean=False,connect=False):
    field=metric+'_mean' if mean else metric
    y=z[field].to_numpy(float)*scale
    err=None
    if metric+'_low' in z:
        err=np.maximum(0,np.vstack([z[field]-z[metric+'_low'],z[metric+'_high']-z[field]])*scale)
    elif metric=='bias':
        err=1.96*z.bias_MCSE.to_numpy(float)
    ax.errorbar(x,y,yerr=err,fmt=MARKERS[model],linestyle='-' if connect else 'none',
                color=PALETTE[model],ms=4.8,mew=.5,markeredgecolor='white',
                capsize=2.3,elinewidth=1.05,linewidth=1.1,zorder=3)


def export(fig,name):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    bounds=fig.bbox
    outside=[]
    inactive_ticks=set()
    for ax in fig.axes:
        for axis in [ax.xaxis,ax.yaxis]:
            lo,hi=sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if not lo-1e-10<=tick.get_loc()<=hi+1e-10:
                    inactive_ticks.update([id(tick.label1),id(tick.label2)])
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip() or id(t) in inactive_ticks:
            continue
        b=t.get_window_extent(renderer)
        if b.x0 < -1 or b.y0 < -1 or b.x1 > bounds.x1+1 or b.y1 > bounds.y1+1:
            outside.append(t.get_text())
    QA.append({'figure':name,'width_mm':round(fig.get_figwidth()*25.4,2),
               'height_mm':round(fig.get_figheight()*25.4,2),'out_of_canvas_text':outside})
    fig.savefig(OUT/(name+'.pdf'),metadata={'Title':name,'Author':'','Creator':'Matplotlib'})
    fig.savefig(OUT/'editable'/(name+'.svg'),metadata={'Title':name,'Creator':'Matplotlib'})
    fig.savefig(OUT/'previews'/(name+'.png'),dpi=300)
    if TIFF:
        from PIL import Image
        (OUT/'tiff').mkdir(exist_ok=True)
        path=OUT/'tiff'/(name+'.tiff')
        fig.savefig(path,dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        with Image.open(path) as source:
            rgb=source.convert('RGB')
        rgb.save(path,compression='tiff_lzw',dpi=(600,600))
        rgb.close()
    plt.close(fig)
    print(name, 'OK' if not outside else 'CHECK '+repr(outside),flush=True)


def fig1():
    z=read('count_summary');z=z[z.de_fraction.gt(0)]
    cases=[(0,0),(0,.5),(.3,0),(.3,.5)]
    fig,axs=plt.subplots(3,4,figsize=(WIDTH,6.1),sharex=True,sharey='row')
    fig.subplots_adjust(left=.10,right=.98,bottom=.10,top=.84,hspace=.35,wspace=.19)
    for c,(sd,mean) in enumerate(cases):
        part=z[z.cage_log2_sd.eq(sd)&z.cage_mean_amplitude.eq(mean)].set_index('model')
        axs[0,c].set_title(f'Cage SD {sd:g}\nMean shift {mean:g}',fontsize=9.4,pad=12)
        for r,metric in enumerate(['FDR','power','coverage']):
            ax=axs[r,c];tidy(ax)
            for j,m in enumerate(COUNT_MODELS):
                points(ax,[j],part.loc[[m]],metric,m,100,True)
            ax.set_xticks([]);ax.set_xlim(-.55,4.55)
            if metric!='power':ref(ax,5 if r==0 else 95)
            if c==0:
                panel(ax,'abc'[r]);ax.set_ylabel(['FDR (%)','Power (%)','Coverage (%)'][r])
    axs[0,0].set_ylim(0,70);axs[1,0].set_ylim(65,102);axs[2,0].set_ylim(70,102)
    legend(fig,COUNT_MODELS,y=.988,ncol=3)
    fig.text(.5,.025,'200 datasets per condition  |  13,312 genes  |  10% nonzero effects',ha='center',fontsize=8.5,color='#526371')
    export(fig,'Fig1')


def fig2():
    z=read('cage_summary');z=z[z.layout.eq('actual_linear')&z.truth.eq(0)]
    fig,axs=plt.subplots(2,3,figsize=(WIDTH,5.55),sharex=True,sharey='col')
    fig.subplots_adjust(left=.105,right=.98,bottom=.13,top=.83,wspace=.36,hspace=.53)
    for r,shift in enumerate([0,.5]):
        for c,metric in enumerate(['bias','rejection','mean_width']):
            ax=axs[r,c];tidy(ax)
            for j,m in enumerate(MODELS):
                a=z[z['shift'].eq(shift)&z.model.eq(m)].sort_values('rho_truth')
                points(ax,a.rho_truth.to_numpy()+(j-2)*.011,a,metric,m,100 if metric=='rejection' else 1)
            ax.set_xticks([0,.2],['0','0.20']);ax.set_xlim(-.065,.265)
            ax.tick_params(labelbottom=True)
            ax.set_ylabel(['Coefficient bias','Null rejection (%)','Mean interval width'][c])
            ax.set_xlabel('Residual correlation')
            panel(ax,chr(97+3*r+c),f'Mean shift {shift:g}')
            if c<2:ref(ax,0 if c==0 else 5)
    legend(fig,MODELS,y=.98,ncol=3)
    export(fig,'Fig2')


def fig3():
    z=read('cage_summary');z=z[z.layout.isin(LAYOUTS)&z.rho_truth.eq(.2)&z['shift'].eq(0)]
    fig,axs=plt.subplots(2,2,figsize=(WIDTH,5.9))
    fig.subplots_adjust(left=.12,right=.98,bottom=.13,top=.83,wspace=.36,hspace=.55)
    for k,(ax,metric,truth,label) in enumerate(zip(axs.flat,['rejection','coverage','rejection','mean_width'],[0,0,.3,.3],['Null rejection (%)','Coverage (%)','Power (%)','Mean interval width'])):
        tidy(ax)
        for j,m in enumerate(MODELS):
            a=z[z.truth.eq(truth)&z.model.eq(m)].set_index('layout').reindex(LAYOUTS)
            points(ax,np.arange(3)+(j-2)*.09,a,metric,m,100 if metric!='mean_width' else 1)
        ax.set_xticks(range(3),['Pure','Half mixed','Fully mixed']);ax.set_xlim(-.4,2.4)
        ax.set_ylabel(label);panel(ax,'abcd'[k])
        if truth==0:ref(ax,5 if metric=='rejection' else 95)
    legend(fig,MODELS,y=.98,ncol=3)
    fig.text(.5,.025,'Pure cages: cage-fixed and Mundlak coefficients are not estimable',ha='center',fontsize=8.4,color='#526371')
    export(fig,'Fig3')


def fig4():
    z=read('general_calibration');z=z[z.status.eq('regular')&z.mean_case.eq('quality_zero')]
    fig,axs=plt.subplots(2,3,figsize=(WIDTH,5.4),sharex=True,sharey='row')
    fig.subplots_adjust(left=.105,right=.98,bottom=.15,top=.82,wspace=.25,hspace=.36)
    for c,law in enumerate(LAWS):
        for r,metric in enumerate(['rejection','coverage']):
            ax=axs[r,c];tidy(ax)
            for m in ['model_based','HC3']:
                a=z[z.error_law.eq(law)&z.method.eq(m)]
                ax.scatter(a.max_leverage,a[metric]*100,s=24,c=PALETTE[m],marker=MARKERS[m],alpha=.8,edgecolors='white',linewidths=.3,zorder=3)
                if law=='heteroskedastic_gaussian':
                    q=a[a.n.eq(80)&a.leverage.eq('ordinary')&a.weights.eq('uniform')&np.isclose(a.partial_correlation,.3)]
                    ax.scatter(q.max_leverage,q[metric]*100,s=78,facecolors='none',edgecolors=INK,linewidths=.8,zorder=4)
            ref(ax,5 if r==0 else 95)
            if c==0:ax.set_ylabel('Null rejection (%)' if r==0 else 'Coverage (%)')
            if r==1:ax.set_xlabel('Maximum leverage')
            panel(ax,chr(97+r*3+c),LAW_LABELS[c] if r==0 else None)
            ax.set_xlim(.05,1.02);ax.set_xticks([.1,.5,.9])
    legend(fig,['model_based','HC3'],y=.99)
    fig.text(.5,.025,'Outlined points: 80 observations, ordinary leverage, equal weights, partial r = 0.3',ha='center',fontsize=8.2)
    export(fig,'Fig4')


def fig5():
    h=read('clinical_pairs').dropna();s=h[h.q_joint_base.lt(.05)].sort_values('q_joint_base')
    fig=plt.figure(figsize=(WIDTH,5.8));gs=fig.add_gridspec(2,2,width_ratios=[1,1],left=.105,right=.98,bottom=.18,top=.89,wspace=.36,hspace=.65)
    a=fig.add_subplot(gs[:,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,1])
    a.hlines(range(40),s.q_joint_base,s.q_joint_quality,color='#C8D5DD',lw=.9)
    a.scatter(s.q_joint_base,range(40),c=BLUE,s=18,marker='o',zorder=3)
    a.scatter(s.q_joint_quality,range(40),c=ORANGE,s=18,marker='s',zorder=3)
    ref(a,.05,True);a.set(xlabel='Joint q value',ylabel='40 base-significant records',yticks=[],xlim=(0,.21),ylim=(-1.5,40.5));a.set_xticks([0,.05,.1,.15,.2]);tidy(a,'x');panel(a,'a','Threshold transitions')
    b.scatter(h.beta_base,h.beta_quality,c='#A8BAC7',s=4,alpha=.45,rasterized=True)
    b.scatter(s.beta_base,s.beta_quality,c=BLUE,s=18,edgecolors='white',linewidths=.3)
    bounds=h[['beta_base','beta_quality']].to_numpy();lo=np.min(bounds);hi=np.max(bounds)
    b.plot([lo,hi],[lo,hi],ls='--',color='#87949F',lw=.8)
    b.set(xlabel='Base coefficient (log2)',ylabel='Adjusted coefficient (log2)');tidy(b);panel(b,'b','All 3,744 clinical records')
    bins=[0,.01,.02,.03,.04,.05,.10,.20];names=['<.01','.01-.02','.02-.03','.03-.04','.04-.05','.05-.10','.10-.20']
    for k,(field,col) in enumerate([('q_joint_base',BLUE),('q_joint_quality',ORANGE)]):
        v=np.histogram(s[field],bins)[0];x=np.arange(7)+(k-.5)*.34
        c.bar(x,v,width=.31,color=col,zorder=2)
        for xx,yy in zip(x,v):
            if yy:c.text(xx,yy+.8,str(yy),ha='center',fontsize=8.2,color=col)
    c.set_xticks(range(7),names,rotation=50,ha='right');c.set(ylabel='Records',ylim=(0,29));tidy(c);panel(c,'c','Selected-record q bins')
    handles=[Line2D([],[],color=BLUE,marker='o',ls='',label='Base'),Line2D([],[],color=ORANGE,marker='s',ls='',label='Quality adjusted')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.025),ncol=2)
    export(fig,'Fig5')


def fig6():
    z=read('family_display');sizes=read('family_sizes');sizes=sizes[sizes.collection.eq('human_17064')]
    scopes=['original','by_contrast','by_context','by_context_contrast']
    names=['Joint family','By contrast','By context','By context + contrast']
    labels=[]
    for sc,n in zip(scopes,names):
        k=sizes[sizes.scope.eq(sc)&sizes.denominator.eq('planned')].planned.unique()
        assert len(k)==1
        labels.append(n+'\n'+f'{int(k[0]):,} positions / family')
    fig,axs=plt.subplots(1,2,figsize=(WIDTH,4.2),sharey=True)
    fig.subplots_adjust(left=.30,right=.97,bottom=.18,top=.83,wspace=.31)
    for j,(branch,col) in enumerate([('base',BLUE),('extended',ORANGE)]):
        a=z[z.record_type.eq('discoveries')&z.denominator.eq('planned')&z.branch.eq(branch)].set_index('scope').reindex(scopes)
        y=np.arange(4)+(j-.5)*.28
        axs[0].barh(y,a.significant,height=.24,color=col)
        for yy,v in zip(y,a.significant):axs[0].text(v+3,yy,f'{v:.0f}',va='center',fontsize=8.3,color=col)
    for j,(den,col) in enumerate([('planned',TEAL),('available_only',PURPLE)]):
        a=z[z.record_type.eq('retention')].set_index('correction')
        v=[a.loc['q__'+s+'__'+den,'adjusted_significant'] for s in scopes]
        y=np.arange(4)+(j-.5)*.28;axs[1].barh(y,v,height=.24,color=col)
        for yy,val in zip(y,v):axs[1].text(val+.7,yy,f'{val:.0f}',va='center',fontsize=8.3,color=col)
    axs[0].set_yticks(range(4),labels);axs[0].invert_yaxis()
    for ax in axs:tidy(ax,'x');ax.set_xlabel('Records')
    axs[0].set_xlim(0,185);axs[1].set_xlim(0,45)
    panel(axs[0],'a','Clinical discoveries');panel(axs[1],'b','Original 40 retained')
    axs[0].legend(handles=[Rectangle((0,0),1,1,color=BLUE,label='Base'),Rectangle((0,0),1,1,color=ORANGE,label='Adjusted')],loc='upper left',bbox_to_anchor=(-.02,-.18),ncol=2)
    axs[1].legend(handles=[Rectangle((0,0),1,1,color=TEAL,label='Planned'),Rectangle((0,0),1,1,color=PURPLE,label='Available')],loc='upper left',bbox_to_anchor=(-.12,-.18),ncol=2)
    export(fig,'Fig6')


def fig7():
    z=read('mouse_summary');ratios=read('mouse_se_ratios');comp=read('mouse_comparisons')
    w=z.pivot(index='spec',columns='model',values='significant_background').reindex(SPECS)
    fig=plt.figure(figsize=(WIDTH,6.1));gs=fig.add_gridspec(3,2,height_ratios=[1.35,1,1],left=.105,right=.98,bottom=.10,top=.84,hspace=.80,wspace=.65)
    a=fig.add_subplot(gs[0,:]);b=fig.add_subplot(gs[1:,0]);c=fig.add_subplot(gs[1,1]);d=fig.add_subplot(gs[2,1])
    pos=b.get_position();b.set_position([.23,pos.y0,.23,pos.height])
    for j,m in enumerate(ARMS):a.bar(np.arange(6)+(j-3)*.105,w[m],width=.095,color=PALETTE[m],zorder=2)
    a.set_xticks(range(6),SPEC_LABELS);a.set_ylabel('BH discoveries');tidy(a);panel(a,'a','Seven matched analysis arms')
    ms=[m for m in ARMS if m!='baseline']
    for j,m in enumerate(ms):
        vals=ratios[ratios.model.eq(m)].set_index('spec').reindex(SPECS).SE_ratio
        b.scatter(vals,j+np.linspace(-.15,.15,len(vals)),c=PALETTE[m],s=30,marker=MARKERS[m],edgecolors='white',linewidths=.4,zorder=3)
    b.set_yticks(range(6),[LABELS[m] for m in ms]);b.invert_yaxis();b.set_xlabel('Median paired SE / no cage');ref(b,1,True);tidy(b,'x');panel(b,'b','Precision by specification')
    for i,s in enumerate(SPECS):c.scatter(w.loc[s,'cage_fixed'],w.loc[s,'mundlak_native'],color=PURPLE,s=35,marker='o',edgecolors='white',linewidths=.4)
    lim=[min(w.cage_fixed.min(),w.mundlak_native.min())*.9,max(w.cage_fixed.max(),w.mundlak_native.max())*1.1]
    c.plot(lim,lim,'--',color='#929EA7',lw=.8);c.set(xlabel='Cage-fixed discoveries',ylabel='Mundlak discoveries');tidy(c);panel(c,'c')
    for j,m in enumerate(['cage_fixed','dream_native']):
        vals=comp[comp.model.eq('mundlak_native')&comp.reference.eq(m)].median_SE_ratio
        d.scatter(vals,j+np.linspace(-.13,.13,len(vals)),color=PALETTE[m],s=30,marker=MARKERS[m],edgecolors='white',linewidths=.4)
    d.set_yticks([0,1],['Cage fixed','dream native']);d.set_xlabel('Median Mundlak / comparator SE');ref(d,1,True);tidy(d,'x');panel(d,'d');d.set_ylim(-.5,1.5)
    legend(fig,ARMS,y=.995,ncol=4)
    export(fig,'Fig7')


def bxp(ax, rows, positions, colors_, horizontal=False, width=.48):
    stats=[dict(med=r['median'],q1=r['q25'],q3=r['q75'],whislo=r.get('whisker_low',r.get('p05')),whishi=r.get('whisker_high',r.get('p95')),fliers=[]) for r in rows]
    q=ax.bxp(stats,positions=positions,widths=width,showfliers=False,patch_artist=True,
             orientation='horizontal' if horizontal else 'vertical',manage_ticks=False,
             medianprops={'color':INK,'linewidth':1},boxprops={'linewidth':.7},
             whiskerprops={'color':'#81909B','linewidth':.7},capprops={'color':'#81909B','linewidth':.7})
    for p,col in zip(q['boxes'],colors_):p.set(facecolor=colors.to_rgba(col,.23),edgecolor=col)


def fig8():
    z=read('precision_components');geom=read('geometry')
    fig,axs=plt.subplots(2,2,figsize=(WIDTH,5.2))
    fig.subplots_adjust(left=.21,right=.98,bottom=.15,top=.88,wspace=.75,hspace=.62)
    factors=[('SE_ratio',BLUE,'SE ratio'),('unscaled_ratio',ORANGE,'Unscaled SE'),('posterior_SD_ratio',TEAL,'Posterior SD')]
    for j,(field,col,label) in enumerate(factors):
        a=z[z.metric.eq(field)].set_index('spec').reindex(SPECS)
        axs[0,0].errorbar(a['median'],np.arange(6)+(j-1)*.19,xerr=[a['median']-a.q25,a.q75-a['median']],fmt=['o','D','s'][j],ms=4.5,color=col,capsize=2,label=label)
    a=z[z.metric.eq('reweight_factor')].set_index('spec').reindex(SPECS)
    bxp(axs[0,1],a.to_dict('records'),range(6),SPEC_COLORS,True)
    for ax in axs[0]:ax.set_yticks(range(6),[SPEC_SHORT[s] for s in SPECS]);ax.invert_yaxis();ref(ax,1,True);tidy(ax,'x')
    axs[0,0].set_xlabel('Median and gene IQR');axs[0,1].set_xlabel('Native / shared fixed-effect SE')
    for ax,metric in zip(axs[1],['genotype_residual_SS','residual_df']):
        for j,m in enumerate(['baseline','cage_fixed']):
            a=geom[geom.branch.eq(m)].set_index('spec').reindex(SPECS[:3])
            ax.barh(np.arange(3)+(j-.5)*.30,a[metric],height=.27,color=PALETTE[m],label=LABELS[m])
        ax.set_yticks(range(3),['Linear age','Age factor','Age + slide']);ax.invert_yaxis();tidy(ax,'x')
    axs[1,0].set_xlabel('Residualized genotype sum of squares');axs[1,1].set_xlabel('Residual degrees of freedom')
    for i,ax in enumerate(axs.flat):panel(ax,'abcd'[i])
    fig.legend(handles=[Line2D([],[],color=c,marker=m,ls='',label=l) for (_,c,l),m in zip(factors,['o','D','s'])],loc='upper center',bbox_to_anchor=(.5,.99),ncol=3)
    handles, labels = axs[1,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center',bbox_to_anchor=(.5,.012),ncol=2)
    export(fig,'Fig8')


def figs1():
    z=read('donors');fig,axs=plt.subplots(1,3,figsize=(WIDTH,3.3))
    fig.subplots_adjust(left=.09,right=.98,bottom=.22,top=.87,wspace=.4)
    for k,(ax,field,label) in enumerate(zip(axs,['PMI','RIN','pH'],['PMI (h)','RIN','CSF pH'])):
        vals=[z[z.clinical_group.eq(g)][field].dropna().to_numpy() for g in ['No_dementia','Dementia']]
        q=ax.boxplot(vals,positions=[0,1],widths=.5,showfliers=False,patch_artist=True,medianprops={'color':INK})
        for p,col in zip(q['boxes'],[BLUE,ORANGE]):p.set(facecolor=colors.to_rgba(col,.16),edgecolor=col)
        for j,(v,col) in enumerate(zip(vals,[BLUE,ORANGE])):
            jitter=((np.arange(len(v))*0.61803398875)%1-.5)*.26
            ax.scatter(j+jitter,v,c=col,s=13,alpha=.7,edgecolors='white',linewidths=.25,zorder=3)
        ax.set_xticks([0,1],['No dementia\nn = 41','Dementia\nn = 39']);ax.set_ylabel(label);tidy(ax);panel(ax,'abc'[k])
    export(fig,'FigS1')


def figs2():
    z=read('astrocyte');contrasts=list(z.contrast.unique());genes=['GFAP','GSDMD','IGFBP5']
    fig,axs=plt.subplots(1,3,figsize=(WIDTH,3.9),sharey=True)
    fig.subplots_adjust(left=.13,right=.98,bottom=.25,top=.76,wspace=.30)
    titles=['Clinical status','Clinical status\nBraak adjusted','Braak stage']
    for c,contrast in enumerate(contrasts):
        ax=axs[c]
        for j,(m,col,marker) in enumerate([('base',BLUE,'o'),('quality',ORANGE,'s')]):
            a=z[z.contrast.eq(contrast)&z.model.eq(m)].set_index('gene_symbol').reindex(genes)
            ax.errorbar(a.beta,np.arange(3)+(j-.5)*.2,xerr=[a.beta-a.low,a.high-a.beta],fmt=marker,ms=5,color=col,elinewidth=1.2,capsize=2)
        ref(ax,0,True);tidy(ax,'x');ax.set_yticks(range(3),genes,fontstyle='italic');ax.set_xlabel('Coefficient (log2)');panel(ax,'abc'[c],titles[c])
    axs[0].invert_yaxis()
    fig.legend(handles=[Line2D([],[],color=BLUE,marker='o',label='Base'),Line2D([],[],color=ORANGE,marker='s',label='Quality adjusted')],loc='lower center',bbox_to_anchor=(.5,.03),ncol=2)
    export(fig,'FigS2')


def figs3():
    z=read('mouse_consistency');f=read('mouse_focal');v=z.significant.value_counts().reindex(range(19),fill_value=0)
    fig,axs=plt.subplots(1,2,figsize=(WIDTH,3.5),gridspec_kw={'width_ratios':[1.4,1]})
    fig.subplots_adjust(left=.09,right=.98,bottom=.22,top=.85,wspace=.35)
    axs[0].bar(v.index,v,color=BLUE,width=.8);axs[0].set(xlabel='Significant specifications (of 18)',ylabel='Genes',xticks=[0,3,6,9,12,15,18]);tidy(axs[0]);panel(axs[0],'a')
    f['discovery']=f.q_three_model_joint.lt(.05);a=f.groupby('gene_symbol').discovery.sum().reindex(['GFAP','GSDMD','IGFBP5','EMP3'])
    axs[1].bar(range(4),a,color=[TEAL,ORANGE,PURPLE,GRAY],width=.65)
    for i,vv in enumerate(a):axs[1].text(i,vv+.5,str(vv),ha='center',fontsize=9)
    axs[1].set_xticks(range(4),['Gfap','Gsdmd','Igfbp5','Emp3'],fontstyle='italic');axs[1].set(ylabel='Significant specifications',ylim=(0,20));tidy(axs[1]);panel(axs[1],'b')
    axs[1].set_yticks([0,5,10,15,20])
    export(fig,'FigS3')


def figs4():
    z=read('icc_density');s=read('icc_strata');bins=['<0.01','0.01-0.10','0.10-0.20','>=0.20']
    bins=list(s.ICC_bin.drop_duplicates())
    # Preserve the source's semantic ICC-bin ordering, independent of row order.
    bins=sorted(bins,key=lambda x:0 if x.startswith('<') else 1 if x.startswith('0.01') else 2 if x.startswith('0.10') else 3)
    fig,axs=plt.subplots(4,3,figsize=(WIDTH,8.5))
    fig.subplots_adjust(left=.10,right=.98,bottom=.15,top=.94,wspace=.40,hspace=.75)
    common_max=0
    hb=[]
    for i,spec in enumerate(SPECS):
        row=0 if i<3 else 2;col=i%3;a=z[z.spec.eq(spec)]
        ax=axs[row,col];h=ax.hexbin(a.ICC,a.SE_ratio,gridsize=28,mincnt=1,cmap='YlGnBu',norm=colors.LogNorm(vmin=1),linewidths=0,rasterized=True);hb.append(h)
        common_max=max(common_max,float(h.get_array().max()));ax.set(xlabel='Conditional ICC',ylabel='Shared dream / no-cage SE');panel(ax,chr(97+i),SPEC_SHORT[spec]);ref(ax,1)
        ax=axs[row+1,col]
        for policy,color,marker in [('shared','#76B9AD','o'),('native',TEAL,'s')]:
            g=s[s.spec.eq(spec)&s.weight_policy.eq(policy)].set_index('ICC_bin').reindex(bins)
            ax.plot(range(4),g.median_SE_ratio,marker=marker,color=color,label=policy,ms=4.5)
        ax.set_xticks(range(4),['<.01','.01-.10','.10-.20','>=.20'],rotation=40,ha='right');ax.set_ylabel('Median SE ratio');ref(ax,1);tidy(ax)
    for h in hb:h.set_norm(colors.LogNorm(vmin=1,vmax=common_max))
    cax=fig.add_axes([.58,.06,.33,.012]);cb=fig.colorbar(hb[0],cax=cax,orientation='horizontal');cb.set_label('Genes per hexagon (log scale)',fontsize=8.2);cb.ax.tick_params(labelsize=8.2)
    fig.legend(handles=[Line2D([],[],color='#76B9AD',marker='o',label='Shared weights'),Line2D([],[],color=TEAL,marker='s',label='Native weights')],loc='lower left',bbox_to_anchor=(.04,.025),ncol=1)
    export(fig,'FigS4')


def figs5():
    z=read('cage_summary');z=z[z.rho_truth.eq(.2)]
    layouts=['actual_linear']+LAYOUTS
    fig,axs=plt.subplots(2,2,figsize=(WIDTH,5.9),sharey='row')
    fig.subplots_adjust(left=.10,right=.98,bottom=.12,top=.82,wspace=.27,hspace=.52)
    for r,truth in enumerate([0,.3]):
        for c,shift in enumerate([0,.5]):
            ax=axs[r,c]
            for j,m in enumerate(MODELS):
                a=z[z.truth.eq(truth)&z['shift'].eq(shift)&z.model.eq(m)].set_index('layout').reindex(layouts)
                points(ax,np.arange(4)+(j-2)*.09,a,'rejection',m,100)
            ax.set_xticks(range(4),['Observed','Pure','Half\nmixed','Fully\nmixed']);ax.set_ylabel('Null rejection (%)' if r==0 else 'Power (%)');ax.set_xlim(-.4,3.4);tidy(ax);panel(ax,chr(97+r*2+c),f'Mean shift {shift:g}')
            if r==0:ref(ax,5)
    legend(fig,MODELS,y=.99,ncol=3)
    export(fig,'FigS5')


def figs6():
    z=read('human_contexts').sort_values('context');y=np.arange(len(z))
    fig,axs=plt.subplots(1,3,figsize=(WIDTH,6.5),sharey=True)
    fig.subplots_adjust(left=.25,right=.98,bottom=.13,top=.88,wspace=.28)
    for ax,field,col,title in zip(axs[:2],['log_abs_beta_ratio','minus_log_SE_ratio'],[BLUE,ORANGE],['Coefficient change','Precision change']):
        m=z[field+'_median'];ax.errorbar(m,y,xerr=[m-z[field+'_q25'],z[field+'_q75']-m],fmt='o',ms=4,color=col,capsize=2);ax.set_xlabel('Log ratio' if col==BLUE else 'Negative log SE ratio');ref(ax,0,True);tidy(ax,'x')
    axs[0].set_yticks(y,z.context.str.replace('_',' ',regex=False));axs[0].invert_yaxis()
    axs[2].barh(y,z.available,color='#BCD4E5',height=.65,label='Available')
    axs[2].barh(y,z.selected,color=PURPLE,height=.65,label='Base significant')
    axs[2].set_xlabel('Records');tidy(axs[2],'x');axs[2].legend(loc='lower left',bbox_to_anchor=(-.3,-.15),fontsize=8.2)
    for i,ax in enumerate(axs):panel(ax,'abc'[i],['Coefficient','Precision','Tested records'][i])
    export(fig,'FigS6')


def figs7():
    d=read('design_display').set_index('spec').reindex(SPECS);z=read('design_errors').set_index('spec').reindex(SPECS)
    fig,axs=plt.subplots(1,2,figsize=(WIDTH,4.2))
    fig.subplots_adjust(left=.11,right=.98,bottom=.21,top=.78,wspace=.65)
    for i,(_,r) in enumerate(d.iterrows()):
        axs[0].errorbar(r.predicted,r.observed_median,yerr=[[r.observed_median-r.observed_q25],[r.observed_q75-r.observed_median]],fmt='o' if i<3 else 's',color=SPEC_COLORS[i],ms=5,capsize=2,label=SPEC_SHORT[SPECS[i]])
    lo=min(d.predicted.min(),d.observed_q25.min())*.95;hi=max(d.predicted.max(),d.observed_q75.max())*1.05
    axs[0].plot([lo,hi],[lo,hi],'--',color='#87949F');axs[0].set(xlabel='Design-only SE ratio',ylabel='Observed SE ratio');tidy(axs[0]);panel(axs[0],'a')
    bxp(axs[1],z.to_dict('records'),range(6),SPEC_COLORS,True);axs[1].set_yticks(range(6),[SPEC_SHORT[s] for s in SPECS]);axs[1].invert_yaxis();ref(axs[1],0,True);axs[1].set_xlabel('Log observed / predicted');tidy(axs[1],'x');panel(axs[1],'b')
    handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.99),ncol=3,fontsize=8.2)
    export(fig,'FigS7')


def figs8():
    z=read('donor_calibration');z=z[z.error_law.eq('matched_gaussian')]
    fig,axs=plt.subplots(1,2,figsize=(WIDTH,3.8))
    fig.subplots_adjust(left=.105,right=.98,bottom=.22,top=.83,wspace=.38)
    a=z[z.mean_case.eq('quality_nonzero_delta_zero')]
    for m in ['model_based','HC3']:
        q=a[a.method.eq(m)];axs[0].scatter(q.joint_F_rejection*100,q.rejection*100,s=17,color=PALETTE[m],marker=MARKERS[m],alpha=.65,edgecolors='white',linewidths=.2)
    ref(axs[0],5);axs[0].set(xlabel='Joint quality F-test rejection (%)',ylabel='Paired-test rejection (%)');tidy(axs[0]);panel(axs[0],'a')
    a=z[z.mean_case.eq('quality_zero')&z.method.eq('HC3')];axs[1].scatter(a.max_leverage,a.RMS_SE/a.empirical_SD,s=18,color=ORANGE,alpha=.7,edgecolors='white',linewidths=.25)
    ref(axs[1],1);axs[1].set(xlabel='Maximum leverage',ylabel='HC3 RMS SE / empirical SD');tidy(axs[1]);panel(axs[1],'b')
    legend(fig,['model_based','HC3'],y=.99)
    export(fig,'FigS8')


def figs9():
    fig,ax=plt.subplots(figsize=(WIDTH,4.5));fig.subplots_adjust(left=.02,right=.98,bottom=.03,top=.95);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    stages=[(.02,.63,.29,.28,BLUE,'1  Match the comparison','Biological units\nTarget coding\nResponse preprocessing'),
            (.355,.63,.29,.28,ORANGE,'2  Compare coefficients','Paired contrast and covariance\nPointwise uncertainty\nTail and family-level error'),
            (.69,.63,.29,.28,TEAL,'3  Locate information','SE components and weights\nWithin-cage overlap\nRank and degrees of freedom')]
    for x,y,w,h,col,title,body in stages:
        ax.add_patch(Rectangle((x,y),w,h,facecolor=colors.to_rgba(col,.07),edgecolor=colors.to_rgba(col,.4),linewidth=.8))
        ax.add_patch(Rectangle((x,y+h-.015),w,.015,color=col))
        ax.text(x+.015,y+h-.06,title,ha='left',va='center',fontsize=8.9,fontweight='bold',color=col)
        ax.text(x+.015,y+.105,body,ha='left',va='center',fontsize=8.9,linespacing=1.8)
    for x0,x1 in [(.31,.355),(.645,.69)]:ax.add_patch(FancyArrowPatch((x0,.77),(x1,.77),arrowstyle='-|>',mutation_scale=12,color='#83939F',lw=1.1))
    ax.add_patch(Rectangle((.12,.10),.76,.32,facecolor='#F2F5F7',edgecolor='#BBC8D1',linewidth=.8))
    ax.text(.5,.35,'Interpret the changed discovery set',ha='center',fontsize=11,fontweight='bold')
    ax.text(.5,.25,'Coefficient change  |  Precision change  |  Testing-family change',ha='center',fontsize=9)
    ax.text(.5,.16,'Report the unit, BH family, availability and relevant calibration',ha='center',fontsize=8.8,color='#526371')
    for x in [.165,.5,.835]:ax.add_patch(FancyArrowPatch((x,.63),(x,.42),arrowstyle='-|>',mutation_scale=12,color='#83939F',lw=1.1))
    export(fig,'FigS9')


def figs10():
    z=read('selected_components').sort_values('log_abs_t_ratio')
    fig,axs=plt.subplots(1,2,figsize=(WIDTH,3.8))
    fig.subplots_adjust(left=.105,right=.98,bottom=.22,top=.85,wspace=.4)
    axs[0].scatter(z.log_abs_beta_ratio,z.minus_log_SE_ratio,s=26,c=BLUE,edgecolors='white',linewidths=.4)
    lim=np.nanmax(np.abs(z[['log_abs_beta_ratio','minus_log_SE_ratio']].to_numpy()))*1.1
    axs[0].plot([-lim,lim],[-lim,lim],'--',color='#A9B5BE',lw=.8);axs[0].plot([-lim,lim],[lim,-lim],'--',color='#A9B5BE',lw=.8)
    ref(axs[0],0);ref(axs[0],0,True);axs[0].set(xlabel='Log absolute coefficient ratio',ylabel='Negative log SE ratio');tidy(axs[0]);panel(axs[0],'a')
    x=np.arange(40)
    axs[1].bar(x-.18,z.log_abs_beta_ratio,width=.36,color=BLUE,label='Coefficient')
    axs[1].bar(x+.18,z.minus_log_SE_ratio,width=.36,color=ORANGE,label='Precision')
    axs[1].plot(x,z.log_abs_t_ratio,color=INK,lw=1.1,label='Sum')
    ref(axs[1],0);axs[1].set(xlabel='Record ordered by log t change',ylabel='Log contribution');tidy(axs[1]);panel(axs[1],'b')
    h,l=axs[1].get_legend_handles_labels();fig.legend(h,l,loc='lower center',bbox_to_anchor=(.5,.025),ncol=3)
    export(fig,'FigS10')


def figs11():
    z=read('donor_calibration');z['SE_ratio']=z.RMS_SE/z.empirical_SD
    laws=['matched_gaussian','RIN_heteroskedastic_gaussian','matched_student_t5']
    fig,axs=plt.subplots(2,2,figsize=(WIDTH,6.0))
    fig.subplots_adjust(left=.11,right=.98,bottom=.13,top=.88,wspace=.35,hspace=.53)
    config=[('rejection','quality_zero','Null rejection (%)',5,100),('coverage','quality_zero','Coverage (%)',95,100),('SE_ratio','quality_zero','RMS SE / empirical SD',1,1),('rejection','delta_two_reference_SD','Two-SD power (%)',None,100)]
    for k,(ax,(metric,case,label,reference,scale)) in enumerate(zip(axs.flat,config)):
        for i,law in enumerate(laws):
            for j,m in enumerate(['model_based','HC3']):
                a=z[z.mean_case.eq(case)&z.error_law.eq(law)&z.method.eq(m)][metric].dropna()*scale
                assert len(a)==144, (case, law, m, len(a))
                pos=i+(j-.5)*.32
                q=a.quantile([0,.25,.5,.75,1]).to_numpy()
                bxp(ax,[dict(zip(['p05','q25','median','q75','p95'],q))],[pos],[PALETTE[m]],width=.26)
                jitter=((np.arange(len(a))*.61803398875)%1-.5)*.17
                ax.scatter(pos+jitter,a,s=7,color=PALETTE[m],marker=MARKERS[m],alpha=.28,linewidths=0,zorder=2)
        ax.set_xticks(range(3),['Matched\nGaussian','RIN-dependent\nGaussian','Matched\nStudent t5']);ax.set_ylabel(label);tidy(ax);panel(ax,'abcd'[k]);ax.set_xlim(-.5,2.5)
        if reference is not None:ref(ax,reference)
    legend(fig,['model_based','HC3'],y=.99)
    export(fig,'FigS11')


def figs12():
    z=read('general_calibration');z=z[z.status.eq('regular')&z.error_law.eq('matched_gaussian')]
    fig,axs=plt.subplots(1,2,figsize=(WIDTH,3.8))
    fig.subplots_adjust(left=.10,right=.98,bottom=.22,top=.86,wspace=.34)
    cancel=z[z.mean_case.eq('quality_nonzero_delta_zero')]
    allv=[]
    for _,g in cancel.groupby('design'):
        a=g.set_index('method');v=np.array([a.loc['model_based','joint_F_rejection'],a.loc['model_based','rejection'],a.loc['HC3','rejection']])*100;allv.append(v)
        axs[0].plot(range(3),v,lw=.75,color='#B8C6D0',alpha=.65,zorder=1)
    a=np.asarray(allv)
    for i,col in enumerate([PURPLE,BLUE,ORANGE]):axs[0].scatter(np.full(len(a),i),a[:,i],c=col,s=17,alpha=.7,zorder=3)
    ref(axs[0],5);axs[0].set_xticks(range(3),['Joint F','Paired\nmodel','Paired\nHC3']);axs[0].set_ylabel('Rejection (%)');tidy(axs[0]);panel(axs[0],'a')
    for i,case in enumerate(['delta_one_reference_SD','delta_two_reference_SD']):
        for j,m in enumerate(['model_based','HC3']):
            a=z[z.mean_case.eq(case)&z.method.eq(m)].rejection*100;pos=3*i+j
            axs[1].scatter(pos+np.linspace(-.13,.13,len(a)),a,color=PALETTE[m],marker=MARKERS[m],s=22,alpha=.75,edgecolors='white',linewidths=.2)
            axs[1].hlines(a.median(),pos-.24,pos+.24,color=INK,lw=1.8)
    axs[1].set_xticks([0,1,3,4],['1 SD\nModel','1 SD\nHC3','2 SD\nModel','2 SD\nHC3']);axs[1].set_ylabel('Pointwise power (%)');tidy(axs[1]);panel(axs[1],'b')
    export(fig,'FigS12')


FUNCTIONS={str(i):globals()['fig'+str(i)] for i in range(1,9)}
FUNCTIONS.update({'S'+str(i):globals()['figs'+str(i)] for i in range(1,13)})
from design_figures import fig9, figs13
from calibration_figures import fig10, figs14, figs15
import sys
FUNCTIONS.update({'9': lambda: fig9(sys.modules[__name__]),
                  'S13': lambda: figs13(sys.modules[__name__])})
FUNCTIONS.update({'10': lambda: fig10(sys.modules[__name__]),
                  'S14': lambda: figs14(sys.modules[__name__]),
                  'S15': lambda: figs15(sys.modules[__name__])})


def main():
    global OUT,TIFF
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'output')
    parser.add_argument('--figures',nargs='+',choices=list(FUNCTIONS),default=list(FUNCTIONS))
    parser.add_argument('--tiff',action='store_true',help='Also export losslessly compressed 600-dpi TIFF files')
    args=parser.parse_args();OUT=args.output;TIFF=args.tiff
    for p in [OUT,OUT/'previews',OUT/'editable']:p.mkdir(parents=True,exist_ok=True)
    style();verify_inputs()
    for name in args.figures:
        FUNCTIONS[name]()
        gc.collect()
    (OUT/'layout_checks.json').write_text(json.dumps(QA,indent=2),encoding='utf-8')
    if any(x['out_of_canvas_text'] for x in QA):
        raise RuntimeError('Inspect labels reported in layout_checks.json')


if __name__=='__main__':
    main()
