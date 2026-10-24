"""Public-input supplementary scientific kernels."""
from pathlib import Path
import csv, gzip, io, json, os, re, shutil, subprocess
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
from pipeline.io import require as need, sha, dump as writej
from .common import csvout, readj, reserved_bh as bh
from scipy.io import mmread
from .common import save, digest
import urllib.parse
from .common import Blocked, utc, need, RDIR
import importlib.metadata
from .annotation import POLICY, age_review
SENSITIVITY_R = (RDIR / 'interactions_template.R').read_text(encoding='utf-8')


COMPARATORS = ['End', 'Ex', 'In', 'Mic', 'Oli', 'Opc']


RISK_REVIEW_POLICY = dict(high_predicted_fraction=.30, expected_rate=.06,
                          histogram_bins=100, scope='descriptive_flags_not_validated_cutoffs',
                          threshold_changes_allowed=False, cell_exclusion_changes_allowed=False)


def make_design(m, model):
    need(model in ('matched_base', 'quality_primary', 'quality_quadratic'), 'Unknown design model')
    fields = ['age'] + (['pmi', 'ph_cerebellum'] if model != 'matched_base' else [])
    need({'unit_id', 'case', 'sex', *fields} <= set(m), 'Missing design metadata fields')
    need(m.unit_id.notna().all() and m.unit_id.is_unique, 'Missing/duplicate design units')
    need(m['case'].isin([0, 1]).all() and set(m['case']) == {0, 1}, 'Design requires binary case 0/1 and both groups')
    need(m.sex.isin(['F', 'M']).all(), 'Unknown/missing sex must not be coded as female')
    for field in fields:
        need(np.isfinite(pd.to_numeric(m[field], errors='coerce').to_numpy(float)).all(),
             'Nonfinite design covariate: ' + field)
    x = pd.DataFrame({'Intercept': np.ones(len(m)), 'case': m.case.to_numpy(float)}, index=m.unit_id)
    for col in ['age'] + (['pmi', 'ph_cerebellum'] if model != 'matched_base' else []):
        a = m[col].to_numpy(float); sd = np.std(a, ddof=1)
        # Constant quality variables add no adjustment and are documented separately.
        if sd > 0:
            x[col] = (a-a.mean())/sd
    if m.sex.nunique() > 1:
        x['male'] = (m.sex == 'M').to_numpy(float)
    if model == 'quality_quadratic' and 'age' in x:
        z = x.age.to_numpy()**2
        x['age_squared'] = z-z.mean()
    need(np.isfinite(x.to_numpy()).all(), 'Nonfinite design')
    return x


def reserved_bh(values):
    x = np.asarray(values, float)
    need(x.ndim == 1 and len(x) > 0, 'BH requires a nonempty one-dimensional family')
    need(not np.isinf(x).any(), 'Infinite p-values are invalid, not missing hypotheses')
    valid = np.isfinite(x)
    need(((x[valid] >= 0) & (x[valid] <= 1)).all(), 'Invalid p-values')
    p = np.where(valid, x, 1.); order = np.argsort(p)
    q = np.empty(len(p)); q[order] = np.minimum(1, np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    q[~valid] = np.nan
    return q


def verify_sensitivity_outputs(d, mode, fixed, expected_arms=None):
    """Reject internally consistent statistics assigned to the wrong endpoint or donor set."""
    z=pd.read_csv(d/'fixed_results.csv')
    need(expected_arms is not None, 'Explicit planned arms required')
    expected=expected_arms
    need(set(z.arm)==set(expected) and len(z)==len(fixed)*len(expected), 'Sensitivity arms/slots missing')
    need(not z.duplicated(['arm','gene_symbol']).any(), 'Duplicate sensitivity result keys')
    meta=pd.read_csv(d/'metadata.csv',dtype={'donor':str})
    need(meta.donor.notna().all(), 'Missing sensitivity donor identity')
    n=meta.donor.nunique()
    from scipy.stats import t
    checks=[]
    for arm, sub in z.groupby('arm',sort=False):
        need(sub.gene_symbol.tolist()==fixed.gene_symbol.tolist(),'Sensitivity endpoint order changed')
        need(sub.testable.notna().all() and sub.testable.isin([True,False]).all(), 'Invalid sensitivity testable flags')
        need(sub.n_donors.eq(n).all(), 'Wrong sensitivity donor count')
        full=pd.read_csv(d/f'all_genes_{arm}.csv.gz')
        need(full.gene_symbol.is_unique, 'Duplicate sensitivity full-gene identities')
        intended=fixed.mapping_status.eq('mapped') & fixed.analysis_symbol.isin(full.gene_symbol)
        need(np.array_equal(sub.testable.to_numpy(bool),intended.to_numpy()), 'Sensitivity endpoint/full-gene mapping mismatch')
        ok=sub.P_Value.notna()
        need(np.array_equal(sub.testable.to_numpy(bool),ok),'Sensitivity testability mismatch')
        numeric=['logFC','CI_low','CI_high','SE_moderated','df_total','P_Value']
        need(sub.loc[~ok,numeric].isna().all().all(), 'Untestable sensitivity endpoint has statistics')
        valid=sub[ok]
        need(np.isfinite(valid[['logFC','CI_low','CI_high','SE_moderated','df_total','P_Value']]).all().all(), 'Nonfinite tested result')
        need((valid.SE_moderated>0).all() and (valid.df_total>0).all(), 'Invalid uncertainty')
        need(np.allclose(valid.P_Value,2*t.sf(abs(valid.logFC/valid.SE_moderated),valid.df_total),rtol=1e-7,atol=1e-12),'Sensitivity p mismatch')
        width=t.ppf(.975,valid.df_total)*valid.SE_moderated
        need(np.allclose(valid.CI_low,valid.logFC-width,rtol=1e-7,atol=1e-10) and
             np.allclose(valid.CI_high,valid.logFC+width,rtol=1e-7,atol=1e-10),'Sensitivity CI mismatch')
        lookup=full.set_index('gene_symbol').reindex(fixed.analysis_symbol)
        for left,right in [('logFC','logFC'),('CI_low','CI.L'),('CI_high','CI.R'),('P_Value','P.Value')]:
            values=lookup[right].to_numpy().copy();values[~intended.to_numpy()]=np.nan
            need(np.allclose(sub[left],values,rtol=1e-10,atol=1e-12,equal_nan=True),
                 'Sensitivity endpoint/full-gene statistic mismatch: '+left)
        checks.append(dict(arm=arm,status='pass',donors=n,tested=int(ok.sum()),
                           checks='identity_mapping;full_gene_link;donor_count;missingness;p;CI'))
    save(d/'output_verification.csv',checks)
    return z


def risk_diagnostics(observed, simulated, threshold, out, sample, annotation=None):
    """Audit the fixed screen without choosing another threshold or changing any flag."""
    need({'barcode','score','predicted_doublet'} <= set(observed), 'Missing observed score fields')
    need(observed.barcode.notna().all() and observed.barcode.is_unique, 'Invalid score barcodes')
    scores = observed.score.to_numpy(float)
    simulated = np.asarray(simulated, float)
    flags = observed.predicted_doublet
    need(len(scores)>0 and len(simulated)>0 and np.isfinite(threshold) and 0<threshold<1,
         'Empty scores or invalid automatic threshold')
    need(np.isfinite(scores).all() and np.isfinite(simulated).all() and
         ((scores>=0)&(scores<=1)).all() and ((simulated>=0)&(simulated<=1)).all(), 'Invalid score range')
    need(flags.notna().all() and flags.isin([True,False]).all() and
         np.array_equal(flags.to_numpy(bool),scores>threshold), 'Score/threshold/flag mismatch')
    observed = observed.copy()
    if annotation is not None:
        ann = annotation.loc[annotation.qc_pass, ['barcode','label','total_umi','detected_genes','mito_fraction']].copy()
        need(ann.barcode.is_unique and set(ann.barcode)==set(observed.barcode), 'QC barcode/score coverage differs')
        observed = observed.drop(columns=[c for c in ann if c!='barcode' and c in observed])
        observed = observed.merge(ann, on='barcode', how='left', validate='one_to_one')
    fraction=float(flags.mean()); detectable=float(np.mean(simulated>threshold))
    issues=[]
    if fraction>RISK_REVIEW_POLICY['high_predicted_fraction']: issues.append('high_predicted_fraction')
    if fraction==0: issues.append('no_observed_scores_above_threshold')
    if detectable in (0.,1.): issues.append('simulated_threshold_outside_interior')
    row=dict(gsm=sample.gsm, fraction=sample.fraction, donor=sample.donor,
             qc_cells=len(scores), predicted_doublets=int(flags.sum()),
             retained_cells=int((~flags.to_numpy(bool)).sum()), threshold=float(threshold),
             predicted_fraction=fraction, simulated_detectable_fraction=detectable,
             observed_over_detectable=fraction/detectable if detectable else None,
             ratio_interpretation='diagnostic_not_certified_doublet_rate',
             algorithm_execution='completed', threshold_review_required=True,
             threshold_validity='unresolved_not_certified',
             diagnostic_flags=';'.join(issues) or 'none_by_fixed_descriptive_rules',
             automatic_threshold_changed=False)
    for population, values in [('observed',scores),('simulated',simulated)]:
        for q in [0,.05,.25,.5,.75,.95,1]: row[f'{population}_q{int(q*100):02d}']=float(np.quantile(values,q))
    bins=np.linspace(0,1,RISK_REVIEW_POLICY['histogram_bins']+1)
    histogram=[]
    for population, values in [('observed',scores),('simulated',simulated)]:
        n,_=np.histogram(values,bins)
        histogram += [dict(population=population,bin_left=bins[i],bin_right=bins[i+1],
                           count=int(v),fraction=float(v/len(values)),threshold=float(threshold)) for i,v in enumerate(n)]
    save(out/'score_histogram.csv',histogram)
    save(out/'observed_scores.csv.gz',observed)
    save(out/'simulated_scores.csv.gz',pd.DataFrame(dict(score=simulated)))
    by_class=[]
    if 'label' in observed:
        for label, group in observed.groupby('label',sort=True):
            by_class.append(dict(label=label,qc_cells=len(group),
                predicted_doublets=int(group.predicted_doublet.sum()),
                predicted_fraction=float(group.predicted_doublet.mean()),
                retained_cells=int((~group.predicted_doublet.to_numpy(bool)).sum())))
        save(out/'class_retention.csv',by_class)
        summaries=[]
        for (label,flag), group in observed.groupby(['label','predicted_doublet'],sort=True):
            for field in ['score','total_umi','detected_genes','mito_fraction']:
                values=group[field].to_numpy(float)
                summaries.append(dict(label=label,predicted_doublet=bool(flag),field=field,n=len(values),
                    median=float(np.median(values)),q25=float(np.quantile(values,.25)),q75=float(np.quantile(values,.75))))
        save(out/'class_QC_score_summary.csv',summaries)
    writej(out/'threshold_diagnostics.json',row)
    return row


def doublet_risk(x, barcodes, keep, out, sample, annotation=None):
    keep = np.asarray(keep)
    need(keep.dtype.kind == 'b' and keep.shape == (len(barcodes),) and x.shape[0] == len(barcodes),
         'Doublet QC mask must be Boolean and aligned with count rows/barcodes')
    ids = np.flatnonzero(keep)
    out.mkdir(parents=True, exist_ok=False)
    if len(ids) < 100:
        row = dict(gsm=sample.gsm, fraction=sample.fraction, donor=sample.donor,
                   status='not_testable_insufficient_qc_cells', qc_cells=len(ids), minimum_qc_cells=100,
                   predicted_doublets=None, retained_cells=None, threshold=None,
                   expected_rate_assumption=.06, seed=290, algorithm_execution='not_run',
                   threshold_review_required=False, threshold_validity='not_available',
                   reason='Too few QC-passing cells for the fixed doublet-risk procedure',
                   boundary='No doublet decision; library excluded from both matched filtering variants')
        writej(out/'specification.json', row)
        # None is unavailable, not an all-singlet or an all-doublet decision.
        return None, row
    import contextlib
    import scrublet
    from threadpoolctl import threadpool_limits
    # Per library only; no artificial mixing of donors, P/V fractions or diagnoses.
    y = x[ids].tocsc()
    y = y[:, np.asarray(y.sum(axis=0)).ravel() > 0].tocsr()
    with (out/'scrublet.log').open('w', encoding='utf-8') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log), threadpool_limits(limits=1):
        sc = scrublet.Scrublet(y, expected_doublet_rate=.06, random_state=290)
        scores, predicted = sc.scrub_doublets(use_approx_neighbors=False, n_prin_comps=30, verbose=True)
    threshold = getattr(sc, 'threshold_', None)
    need(threshold is not None and np.isfinite(threshold) and predicted is not None,
         'Scrublet did not find a usable automatic threshold; no forced threshold')
    scores = np.asarray(scores, float); predicted = np.asarray(predicted)
    need(scores.shape == (len(ids),) and np.isfinite(scores).all() and predicted.dtype.kind == 'b', 'Invalid Scrublet output')
    need(np.array_equal(predicted, scores > threshold), 'Doublet flags do not follow the recorded threshold')
    simulated = np.asarray(sc.doublet_scores_sim_, float)
    need(np.isfinite(simulated).all(), 'Nonfinite simulated doublet scores')
    diagnostic = risk_diagnostics(pd.DataFrame(dict(barcode=np.asarray(barcodes)[ids], score=scores,
        predicted_doublet=predicted)), simulated, threshold, out, sample, annotation)
    result = np.zeros(len(barcodes), bool); result[ids] = ~predicted
    row = dict(gsm=sample.gsm, fraction=sample.fraction, donor=sample.donor,
               status='completed_risk_screen_not_doublet_truth', qc_cells=len(ids),
               predicted_doublets=int(predicted.sum()), retained_cells=int(result.sum()),
               threshold=float(threshold), expected_rate_assumption=.06, seed=290,
               package_version=importlib.metadata.version('scrublet'),
               boundary='automatic_threshold_requires_distribution_review;not_certified_cell_identity',
               **{k:diagnostic[k] for k in ['predicted_fraction','simulated_detectable_fraction',
                  'algorithm_execution','threshold_review_required','threshold_validity','diagnostic_flags']})
    writej(out/'specification.json', row)
    return result, row


def matching_gate(m):
    need(m.donor.notna().all() and m.donor.is_unique and m.unit_id.is_unique, 'Missing/duplicate analysis identity')
    need(m.case.isin([0,1]).all() and min(m.case.eq(0).sum(),m.case.eq(1).sum()) >= POLICY['min_group'], 'Insufficient group coverage')
    need(m[['age','pmi','ph_cerebellum']].notna().all().all(), 'Missing required quality covariate')
    age = age_review(m)
    need(age['pass_gate'], 'Age coverage gate failed: '+str(age))
    X = make_design(m,'quality_primary')
    need(np.linalg.matrix_rank(X) == X.shape[1] and len(X)-X.shape[1] >= POLICY['min_df'], 'Non-estimable donor design')
    return X


def boundary_engine(mode, comparator=None):
    script = SENSITIVITY_R
    if mode in ('interaction','interaction_second_rho'):
        need(comparator in COMPARATORS, 'Comparator outside the fixed six classes')
        # Replace only the explicit class literals/column labels, never an R API name.
        replacements={"'Oli'":repr(comparator), "'Oli_'":repr(comparator+'_'),
            "'Oli_case'":repr(comparator+'_case'), 'Oli=ko':comparator+'=ko',
            "'Ast_minus_Oli_interaction'":repr('Ast_minus_'+comparator+'_interaction'),
            "'AD_minus_CT_Ast_minus_AD_minus_CT_Oli'":repr('AD_minus_CT_Ast_minus_AD_minus_CT_'+comparator)}
        for old,new in replacements.items():
            need(old in script,'R template anchor missing: '+old)
            script=script.replace(old,new)
        script=script.replace('Requires exactly Ast and Oli labels','Requires exactly Ast and the specified comparator labels')
        script=script.replace('Requires paired Ast and Oli','Requires paired Ast and the specified comparator')
        anchor="dc<-duplicateCorrelation(v0,X,block=m$donor);rho<-dc$consensus.correlation"
        need(script.count(anchor)==1,'Correlation template anchor changed')
        script=script.replace(anchor,anchor+'\n rho_first<-rho;rho_second<-NA_real_')
        script=script.replace('data.frame(rho=rho,donors=',
            'data.frame(rho=rho,rho_first=rho_first,rho_second=rho_second,donors=')
        arm='Ast_minus_'+comparator+'_interaction'
        if mode=='interaction_second_rho':
            anchor="v<-voom(dy,X,block=m$donor,correlation=rho,plot=FALSE,span=.5,adaptive.span=TRUE)"
            need(script.count(anchor)==1,'Second-rho voom anchor changed')
            script=script.replace(anchor,anchor+'''
 dc2<-duplicateCorrelation(v,X,block=m$donor);rho_second<-dc2$consensus.correlation
 need(is.finite(rho_second)&&abs(rho_second)<.99,'Unusable second within-donor correlation; no first-rho fallback')
 rho<-rho_second
 v<-voom(dy,X,block=m$donor,correlation=rho,plot=FALSE,span=.5,adaptive.span=TRUE)''')
            script=script.replace(repr(arm),repr(arm+'_second_rho'))
            arm+='_second_rho'
        expected=[arm]
    else:
        raise ValueError('Only the planned interaction sensitivity is supported')
    need('dc$singlet.correlation' not in script and 'dc$consensus.correlation' in script,
         'An R API field was changed by template generation')
    return script, expected
