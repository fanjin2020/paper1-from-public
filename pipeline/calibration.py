"""在SEA-AD实际供体设计上检验配对系数差的条件频率学性质。

模拟的是归一化表达尺度的线性响应，不读取原始计数，也不重新估计voom权重。
检验单位为一个预先给定的目标系数差；不把逐项拒绝率解释为BH后的FDR。
"""
import os
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(name, '1')
from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import pandas as pd
from statistics import NormalDist


def zscore(x):
    """中心化并缩放设计列；常数列返回零，后续秩检查仍会发现退化。"""
    x = np.asarray(x, float)
    s = x.std()
    return (x - x.mean()) / s if s > 0 else np.zeros_like(x)


def wilson(k, n):
    """用固定计划数或明确的有效次数计算Wilson区间。"""
    if n == 0:
        return np.nan, np.nan
    z = NormalDist().inv_cdf(.975)
    p = k / n
    center = (p + z*z/(2*n))/(1+z*z/n)
    half = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return center-half, center+half


def geometry(x0, x1, w, target=1):
    """构造同一响应的两个WLS映射及差值，保留共享响应协方差。

    d是完整模型与基础模型目标行之差。HC3使用完整模型残差和其杠杆值。
    不相减两套独立标准误，因为两个系数来自同一组观测。
    """
    x0, x1, w = np.asarray(x0,float), np.asarray(x1,float), np.asarray(w,float)
    n = len(w)
    if not (len(x0) == len(x1) == n and np.all(np.isfinite(w)) and np.all(w > 0)):
        raise ValueError('设计维度或权重不合法')
    for x in (x0, x1):
        if np.linalg.matrix_rank(x) != x.shape[1]:
            raise ValueError('设计矩阵秩不足')
    if np.linalg.matrix_rank(np.column_stack((x1,x0))) != x1.shape[1]:
        raise ValueError('两个模型不满足嵌套关系')
    a0 = np.linalg.solve(x0.T @ (w[:,None]*x0), x0.T*w)
    a1 = np.linalg.solve(x1.T @ (w[:,None]*x1), x1.T*w)
    d = a1[target]-a0[target]
    unit = np.sum(d*d/w)
    h = np.sum(x1*a1.T, axis=1)
    df = n-x1.shape[1]
    if df <= 0 or h.max() >= 1-1e-10:
        raise ValueError('残差自由度或杠杆值不可用')
    return dict(x0=x0,x1=x1,w=w,a0=a0,a1=a1,d=d,unit=unit,h=h,df=df)


def evaluate(y, g):
    """批量计算配对差值、两类标准误和新增质量项的联合F统计量。"""
    e1 = y-g['x1']@(g['a1']@y)
    e0 = y-g['x0']@(g['a0']@y)
    rss1 = np.sum(g['w'][:,None]*e1**2,axis=0)
    rss0 = np.sum(g['w'][:,None]*e0**2,axis=0)
    delta = g['d']@y
    se_model = np.sqrt(rss1/g['df']*g['unit'])
    se_hc3 = np.sqrt((g['d']**2/(1-g['h'])**2)@(e1**2))
    added = g['x1'].shape[1]-g['x0'].shape[1]
    f = np.maximum(0,(rss0-rss1)/added)/(rss1/g['df'])
    return delta,se_model,se_hc3,f


def build_designs(metadata):
    """按新下载数据中的类群和三项对比建立设计。

    age及病理列的缩放不改变协变量空间。目标列保留原始临床0/1或Braak级数。
    质量列缩放只用于构造可解释的模拟均值，不改变嵌套模型的投影。
    """
    required = ['age','sex','pathology','PMI','RIN','pH']
    for context, frame in metadata.groupby('context',sort=True):
        good = (frame.n_cells >= 20) & frame.clinical_group.isin(['Dementia','No_dementia'])
        d = frame.loc[good].dropna(subset=required).copy()
        if d.donor_id.duplicated().any():
            raise ValueError('同一类群内出现重复供体')
        for contrast in ('clinical_status','clinical_status_adjusted_for_braak','braak_same_clinical_population'):
            target = d.pathology.to_numpy(float) if contrast.startswith('braak') else (d.clinical_group=='Dementia').to_numpy(float)
            cols = [np.ones(len(d)),target]
            if d.age.nunique()>1:
                cols.append(zscore(d.age))
            if d.sex.nunique()>1:
                cols.append((d.sex=='male').to_numpy(float))
            if contrast=='clinical_status_adjusted_for_braak':
                cols.append(zscore(d.pathology))
            x0 = np.column_stack(cols)
            z = np.column_stack([zscore(d[c]) for c in ('PMI','RIN','pH')])
            yield context,contrast,d,x0,np.column_stack((x0,z))


def mean_profiles(g):
    """构造零质量效应、抵消型零差异和非零差异三个已知均值。

    抵消情形的三个质量系数不全为零，但其对目标系数差的投影为零。
    非零差异设为匹配误差模型下的两倍解析SD；不是生物学最小效应。
    """
    z = g['x1'][:,-3:]
    c = g['d']@z
    cn = np.linalg.norm(c)
    if cn < 1e-10 or g['unit'] < 1e-14:
        raise ValueError('差值退化，不能构造局部备择')
    candidate = np.eye(3)[np.argmin(np.abs(c))]
    q = candidate-c*np.dot(candidate,c)/(cn*cn)
    v = z@q
    v = v / np.sqrt(np.mean(v*v))*.5
    shift = z@c * (2*np.sqrt(g['unit'])/(cn*cn))
    # 基础均值完整地落在两个模型共有的列空间，数值不会造成目标差异。
    common = .2*g['x0'][:,1]
    if g['x0'].shape[1]>2:
        common = common+.1*g['x0'][:,2]
    return [('quality_zero',common),('quality_nonzero_delta_zero',common+v),('delta_two_reference_SD',common+shift)]


def run(metadata_path, output, repeats=2500, seed=20261005):
    """执行完整事先写定的网格；保存逐情境结果、几何诊断及运行配置。"""
    from scipy import stats
    output = Path(output)
    output.mkdir(parents=True,exist_ok=True)
    # The core stage receipt controls continuation. An unreceipted summary is
    # incomplete and must be recalculated, never accepted as completed work.
    meta = pd.read_csv(metadata_path)
    config = dict(repeats=repeats,seed=seed,contexts=int(meta.context.nunique()),contrasts=3,
                  weights=['uniform','fixed_nucleus_count_gradient'],
                  error_laws=['matched_gaussian','RIN_heteroskedastic_gaussian','matched_student_t5'],
                  means=['quality_zero','quality_nonzero_delta_zero','delta_two_reference_SD'],
                  interval_level=.95,alpha=.05,testing='pointwise; no BH simulation',
                  weights_are_voom=False,expression_pipeline_reestimated=False,
                  metadata_sha256=hashlib.sha256(Path(metadata_path).read_bytes()).hexdigest())
    (output/'calibration_plan.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
    records,geom, unavailable = [],[],[]
    index = 0
    for design_id,(context,contrast,md,x0,x1) in enumerate(build_designs(meta)):
        for wi,weight_name in enumerate(config['weights']):
            w = np.ones(len(md)) if wi==0 else np.exp(.5*np.clip(zscore(np.log1p(md.n_cells)), -2,2))
            w /= w.mean()
            try:
                if len(md)<10 or md.clinical_group.value_counts().reindex(['Dementia','No_dementia'],fill_value=0).min()<5:
                    raise ValueError('完整病例每组不足5')
                g = geometry(x0,x1,w)
                if g['df']<5:
                    raise ValueError('残差自由度不足5')
                profiles = mean_profiles(g)
            except (ValueError, np.linalg.LinAlgError) as exc:
                unavailable.append(dict(context=context,contrast=contrast,weight=weight_name,reason=str(exc)))
                continue
            geom.append(dict(context=context,contrast=contrast,weight=weight_name,n=len(md),p0=x0.shape[1],p1=x1.shape[1],
                             df=g['df'],max_leverage=g['h'].max(),unit_variance=g['unit'],
                             weight_min=w.min(),weight_max=w.max(),condition=np.linalg.cond(x1*np.sqrt(w[:,None]))))
            for ei,error in enumerate(config['error_laws']):
                rng = np.random.default_rng(np.random.SeedSequence([seed,design_id,wi,ei]))
                var = 1/w
                if ei==1:
                    multiplier = np.exp(.7*np.clip(zscore(md.RIN),-2,2))
                    var = var*multiplier/np.mean(multiplier)
                noise = rng.standard_t(5,size=(len(md),repeats))*np.sqrt(3/5) if ei==2 else rng.standard_normal((len(md),repeats))
                noise *= np.sqrt(var[:,None])
                oracle_sd = np.sqrt(np.sum(g['d']**2*var))
                for mean_name,mu in profiles:
                    truth = float(g['d']@mu)
                    if mean_name!='delta_two_reference_SD':
                        assert abs(truth)<1e-10
                        truth=0.
                    delta,se0,se1,f = evaluate(mu[:,None]+noise,g)
                    f_p = stats.f.sf(f,x1.shape[1]-x0.shape[1],g['df'])
                    for method,se in [('model_based',se0),('HC3',se1)]:
                        valid = np.isfinite(delta)&np.isfinite(se)&(se>1e-12)
                        err = delta[valid]-truth
                        crit=stats.t.ppf(.975,g['df'])
                        covered=np.abs(err)<=crit*se[valid]
                        rejected=np.abs(delta[valid]/se[valid])>crit
                        na=int(valid.sum())
                        cr,rr=int(covered.sum()),int(rejected.sum())
                        clo,chi=wilson(cr,na);rlo,rhi=wilson(rr,na)
                        row=dict(context=context,contrast=contrast,weight=weight_name,error_law=error,mean_case=mean_name,method=method,
                                 planned=repeats,available=na,failures=repeats-na,failure_rate=(repeats-na)/repeats,
                                 true_delta=truth,oracle_SD=oracle_sd,bias=np.mean(err),bias_MCSE=np.std(err,ddof=1)/np.sqrt(na),
                                 empirical_SD=np.std(delta[valid],ddof=1),mean_SE=np.mean(se[valid]),RMS_SE=np.sqrt(np.mean(se[valid]**2)),
                                 coverage=cr/na,coverage_low=clo,coverage_high=chi,coverage_all_planned=cr/repeats,
                                 mean_width=np.mean(2*crit*se[valid]),rejection=rr/na,rejection_low=rlo,rejection_high=rhi,
                                 rejected_all_planned=rr/repeats,rejection_MCSE=np.sqrt(rr/na*(1-rr/na)/na),
                                 joint_F_rejection=float(np.mean(f_p<.05)),max_leverage=g['h'].max(),n=len(md),df=g['df'])
                        records.append(row)
                    index+=1
        print(f'校准设计 {design_id+1}/{config["contexts"]*3} {context} {contrast}',flush=True)
    frame=pd.DataFrame(records)
    if frame.empty:
        raise ValueError('没有可校准设计')
    assert frame[['context','contrast','weight','error_law','mean_case','method']].duplicated().sum()==0
    frame.to_csv(output/'calibration_summary.csv',index=False)
    pd.DataFrame(geom).to_csv(output/'calibration_geometry.csv',index=False)
    pd.DataFrame(unavailable,columns=['context','contrast','weight','reason']).to_csv(output/'calibration_unavailable.csv',index=False)
    info=dict(status='completed',scenarios=index,replicates_per_scenario=repeats,model_evaluations=len(frame)*repeats,
              failures=int(frame.failures.sum()),python=platform.python_version(),numpy=np.__version__,
              conditional_fixed_design=True,clinical_or_biological_validation=False)
    (output/'calibration_status.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
    print(json.dumps(info,indent=2),flush=True)
    return frame


