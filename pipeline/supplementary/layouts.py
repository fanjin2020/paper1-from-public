"""Cage-layout sensitivity, using current public sample metadata."""
import itertools, math
import numpy as np
import pandas as pd
from scipy import linalg, stats
from .common import need
MODELS = ('no_cage', 'cage_fixed', 'cage_correlation')
def synthetic_units(size,mixed_fraction,mutant_fraction,seed):
    n=72;k=n//size;m=round(k*mixed_fraction);target=round(n*mutant_fraction)
    rng=np.random.default_rng(seed)
    possibilities=[pure for pure in range(k-m+1) if m<=target-pure*size<=m*(size-1)]
    if not possibilities:return None
    pure=possibilities[len(possibilities)//2]
    counts=[1]*m;left=target-pure*size-m
    while left:
        available=[i for i,v in enumerate(counts) if v<size-1]
        smallest=min(counts[i] for i in available)
        available=[i for i in available if counts[i]==smallest]
        j=int(rng.choice(available));counts[j]+=1;left-=1
    counts+=([size]*pure+[0]*(k-m-pure));rng.shuffle(counts)
    rows=[]
    for c,num in enumerate(counts):
        g=[1]*num+[0]*(size-num);rng.shuffle(g)
        rows.extend(dict(cage=c,genotype=v) for v in g)
    u=pd.DataFrame(rows)
    need(len(u)==n and u.genotype.sum()==target,'Synthetic layout genotype counts differ')
    need(u.groupby('cage').genotype.nunique().eq(2).sum()==m,'Mixed-cage count differs')
    return u

def designs(u,kind):
    base=pd.DataFrame({'Intercept':np.ones(len(u)),'genotype':u.genotype.to_numpy(float)})
    if kind=='linear':base['age_centered']=u.age-u.age.mean();base['male']=u.male
    if kind in ['age_factor','age_factor_slide']:
        base=pd.concat([base,pd.get_dummies(u.age,prefix='age',drop_first=True,dtype=float)],axis=1)
        base['male']=u.male
    if kind=='age_factor_slide':
        base=pd.concat([base,pd.get_dummies(u.slide,prefix='slide',drop_first=True,dtype=float)],axis=1)
    need(np.linalg.matrix_rank(base)==base.shape[1],'Base design is not full rank')
    nuisance=pd.concat([base.drop(columns='genotype'),pd.get_dummies(u.cage,prefix='cage',drop_first=True,dtype=float)],axis=1)
    q,r,piv=linalg.qr(nuisance.to_numpy(float),mode='economic',pivoting=True)
    rank=np.linalg.matrix_rank(nuisance)
    # A basis removes only redundant nuisance columns, never the genotype contrast.
    basis=nuisance.iloc[:,piv[:rank]]
    g=base.genotype.to_numpy();res=g-q[:,:rank]@(q[:,:rank].T@g)
    ss=float(res@res)
    fixed=pd.concat([base[['genotype']],basis],axis=1)
    n0=base.drop(columns='genotype').to_numpy();r0=g-n0@np.linalg.lstsq(n0,g,rcond=None)[0]
    ss0=float(r0@r0)
    estimable=ss>1e-9 and len(u)-fixed.shape[1]>1
    cage_mean=u.groupby('cage').genotype.transform('mean').to_numpy()
    shift=cage_mean-g.mean()
    info=dict(n_animals=len(u),n_cages=u.cage.nunique(),mutant_fraction=float(g.mean()),
       actual_mixed_fraction=float(u.groupby('cage').genotype.nunique().eq(2).mean()),
       base_columns=base.shape[1],fixed_nuisance_rank=int(rank),
       fixed_redundant_nuisance_columns=nuisance.shape[1]-rank,
       base_df=len(u)-base.shape[1],fixed_df=len(u)-np.linalg.matrix_rank(fixed),
       genotype_SS_base=ss0,genotype_SS_fixed=ss,information_fraction=ss/ss0,
       fixed_estimable=bool(estimable),unit_shift_SD=float(np.std(shift,ddof=1)),
       analytic_no_cage_bias_per_shift=float(np.linalg.lstsq(base,shift,rcond=None)[0][1]))
    return base,fixed,info

def build_layouts(seed, actual):
    layouts={};registry=[]
    actual=actual[['cage','genotype','age','male','slide']].copy().reset_index(drop=True)
    need(not actual.isna().any().any() and actual.genotype.isin([0,1]).all(), 'Invalid current animal metadata')
    need(actual.cage.nunique()>1 and actual.genotype.nunique()==2, 'Missing cage or genotype coverage')
    for kind in ['linear','age_factor','age_factor_slide']:
        name='actual_'+kind;b,f,info=designs(actual,kind)
        layouts[name]=(actual.copy(),b,f);registry.append(dict(layout_id=name,kind='actual',feasible=True,**info))
    for idx,(size,mix,frac) in enumerate(itertools.product([2,4,6],[0.,.5,1.],[1/3,.5])):
        name=f'n72_s{size}_m{int(mix*100):03d}_g{round(frac*72):02d}'
        u=synthetic_units(size,mix,frac,seed+idx)
        row=dict(layout_id=name,kind='synthetic',cage_size=size,requested_mixed_fraction=mix,
                 requested_mutant_fraction=frac,feasible=u is not None)
        if u is None:
            row['reason']='Exact genotype balance and mixed-cage count cannot both be achieved'
        else:
            b,f,info=designs(u,'none');layouts[name]=(u,b,f);row.update(info)
        registry.append(row)
    return layouts,pd.DataFrame(registry)

def wilson(k,n):
    if n==0:return np.nan,np.nan
    p=k/n;z=stats.norm.ppf(.975);d=1+z*z/n
    c=(p+z*z/(2*n))/d;w=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return c-w,c+w

def summarize(raw,scenarios,registry):
    data=raw.merge(scenarios,on='scenario_id',validate='many_to_one')
    rows=[];paired=[]
    for (sid,model),z in data.groupby(['scenario_id','model'],sort=False):
        ok=z.status.eq('ok');a=z.loc[ok];n=len(a);planned=len(z)
        base={k:z.iloc[0][k] for k in ['scenario_id','layout_id','rho','shift','truth']}
        row=dict(base,model=model,n_planned=planned,n_success=n,n_not_estimable=int(z.status.eq('not_estimable').sum()),
                 n_failed=int(z.status.eq('failed').sum()),failure_rate=1-n/planned)
        row['failure_low'],row['failure_high']=wilson(planned-n,planned)
        metrics=['bias','bias_MCSE','empirical_SD','mean_SE','RMS_SE','SE_to_empirical_SD','coverage',
                 'coverage_MCSE','coverage_low','coverage_high','mean_CI_width','rejection_rate',
                 'rejection_MCSE','rejection_low','rejection_high','null_rejection_rate','mean_rho_estimate',
                 'mean_gene_clip_fraction','near_pd_boundary_rate','warning_rate']
        row.update({k:np.nan for k in metrics})
        if n:
            errors=a.beta-a.truth;covered=a.ci_low.le(a.truth)&a.ci_high.ge(a.truth);rej=a.p.lt(.05)
            sd=float(errors.std(ddof=1)) if n>1 else np.nan
            rms=float(np.sqrt(np.mean(a.se**2)))
            row.update(bias=float(errors.mean()),bias_MCSE=sd/math.sqrt(n),empirical_SD=sd,
              mean_SE=float(a.se.mean()),RMS_SE=rms,SE_to_empirical_SD=rms/sd if sd>0 else np.nan,
              coverage=float(covered.mean()),coverage_MCSE=math.sqrt(covered.mean()*(1-covered.mean())/n),
              mean_CI_width=float((a.ci_high-a.ci_low).mean()),rejection_rate=float(rej.mean()),
              rejection_MCSE=math.sqrt(rej.mean()*(1-rej.mean())/n),
              null_rejection_rate=float(rej.mean()) if z.truth.iloc[0]==0 else np.nan)
            row['coverage_low'],row['coverage_high']=wilson(int(covered.sum()),n)
            row['rejection_low'],row['rejection_high']=wilson(int(rej.sum()),n)
        if model=='cage_correlation':
            row['n_rho_finite']=int(np.isfinite(z.rho_estimate).sum())
            row['mean_rho_estimate']=z.rho_estimate.mean()
            rho_values=z.rho_estimate[np.isfinite(z.rho_estimate)]
            row['rho_SD']=rho_values.std(ddof=1)
            row['rho_q025']=rho_values.quantile(.025)
            row['rho_q975']=rho_values.quantile(.975)
            row['negative_rho_fraction']=rho_values.lt(0).mean() if len(rho_values) else np.nan
            row['mean_gene_clip_fraction']=z.gene_correlation_near_clip_fraction.mean()
            values=z.consensus_near_pd_boundary.dropna().astype(str).str.lower().eq('true')
            row['near_pd_boundary_rate']=values.mean() if len(values) else np.nan
            row['warning_rate']=z.warning.fillna('').ne('').mean()
        rows.append(row)
    for sid,z in data.groupby('scenario_id',sort=False):
        anchor=z[z.model.eq('no_cage')]
        for model in MODELS[1:]:
            p=anchor.merge(z[z.model.eq(model)],on='replicate',suffixes=('_a','_b'),validate='one_to_one')
            p=p[p.status_a.eq('ok')&p.status_b.eq('ok')]
            for metric in ['beta','se','ci_width']:
                diff=(p.ci_high_b-p.ci_low_b)-(p.ci_high_a-p.ci_low_a) if metric=='ci_width' else p[metric+'_b']-p[metric+'_a']
                se=diff.std(ddof=1)/math.sqrt(len(diff)) if len(diff)>1 else np.nan
                paired.append(dict(scenario_id=sid,model=model,reference='no_cage',metric=metric,n_paired=len(p),
                  mean_difference=diff.mean(),paired_MCSE=se,low=diff.mean()-1.96*se,high=diff.mean()+1.96*se))
    summary=pd.DataFrame(rows).merge(registry,on='layout_id',validate='many_to_one')
    return summary,pd.DataFrame(paired)
