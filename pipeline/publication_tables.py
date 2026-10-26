"""Publication table definitions derived solely from a completed experiment.

No reference workbook, historical results, expression fitting, or simulation is
used here. Native two-factor and shared-weight three-factor identities are kept
separate; ICC is paired with dream/base, never with fixed/base.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .io import require
from .supplementary.common import reserved_bh

VERSION = 'publication_exports_20261026_v3'
BLOCK_COUNTS = (3,3,11,4,5,6,4,3,4,2,6,11,4,5,4,4,7,3,6,17)
EXPECTED_BLOCKS = {f'S{s}.{i}' for s,n in enumerate(BLOCK_COUNTS,1) for i in range(1,n+1)}
INPUT_DEFINITIONS_URL = 'https://github.com/fanjin2020/paper1-from-public/blob/v2.1.0/docs/DATA.md'
KEY = ['spec', 'gene_symbol']
ICC_EDGES = [-np.inf, .01, .10, .20, np.inf]
ICC_LABELS = ['<0.01', '0.01-0.10', '0.10-0.20', '>=0.20']


def unique(frame, keys, label):
    require(not frame.duplicated(keys).any(), 'Duplicate keys in ' + label)
    return frame


def native_precision(raw, precision):
    a = unique(raw[raw.branch.eq('baseline')], KEY, 'base decomposition').set_index(KEY)
    b = unique(raw[raw.branch.eq('cage_fixed')], KEY, 'fixed decomposition').set_index(KEY)
    require(a.index.difference(b.index).empty and b.index.difference(a.index).empty,
            'Native decomposition gene sets differ')
    b = b.reindex(a.index)
    z = pd.DataFrame(index=a.index)
    z['beta_base'] = a.beta_direct
    z['beta_fixed'] = b.beta_direct
    z['SE_base'] = a.SE_direct_ordinary
    z['SE_fixed'] = b.SE_direct_ordinary
    z['unscaled_ratio'] = b.unscaled_direct / a.unscaled_direct
    z['posterior_variance_ratio'] = (z.SE_fixed / z.SE_base / z.unscaled_ratio)**2
    p = unique(precision, KEY, 'shared decomposition').set_index(KEY)
    require(z.index.difference(p.index).empty and p.index.difference(z.index).empty
            and p.SE_ratio.notna().all(),
            'Missing shared precision records')
    z = z.join(p)
    require(np.allclose(z.SE_fixed/z.SE_base, z.SE_ratio, rtol=1e-8, atol=1e-10),
            'Native fitted SE ratios disagree with recorded decomposition')
    z['SE_shared_weights'] = z.SE_fixed / z.reweight_factor
    z['design_same_weights_unscaled_ratio'] = z.design_factor
    z['shared_weights_posterior_variance_ratio'] = z.posterior_factor**2
    z['reweighting_SE_ratio'] = z.reweight_factor
    z['posterior_SD_ratio'] = np.sqrt(z.posterior_variance_ratio)
    z['observed_ratio'] = z.SE_ratio
    z['identity_ratio'] = z.unscaled_ratio * z.posterior_SD_ratio
    z['common_weight_identity_ratio'] = z.design_factor*z.posterior_factor*z.reweight_factor
    z['design_only_ratio'] = z.unweighted_design_prediction
    z['design_only_log_error'] = np.log(z.SE_ratio/z.design_only_ratio)
    z['weighted_information_log_error'] = np.log(z.SE_ratio/z.design_factor)
    require(np.allclose(z.common_weight_identity_ratio, z.SE_ratio, rtol=1e-8, atol=1e-10),
            'Shared-weight three-factor identity failed')
    return z.reset_index()


def icc_precision(variance, background):
    base = unique(background[background.model.eq('baseline')], KEY, 'background base')
    base = base[KEY+['SE']].rename(columns={'SE': 'SE_base'})
    dr = background[background.model.isin(['dream_shared', 'dream_native'])].copy()
    dr['weight_policy'] = dr.model.str.replace('dream_', '', regex=False)
    dr = dr[KEY+['weight_policy', 'SE', 'beta', 'status']].rename(
        columns={'beta': 'beta_dream', 'status': 'dream_status'})
    dr = unique(dr, KEY+['weight_policy'], 'dream fits').merge(base, on=KEY, validate='many_to_one')
    v = unique(variance, KEY+['weight_policy'], 'variance fits')
    require(len(v) == len(dr), 'ICC/dream record counts differ')
    z = v.merge(dr, on=KEY+['weight_policy'], how='outer', validate='one_to_one', indicator=True)
    require(z['_merge'].eq('both').all(), 'ICC/dream gene sets differ')
    z = z.drop(columns='_merge')
    z['SE_ratio'] = z.SE/z.SE_base
    z['SE_ratio_definition'] = 'dream_' + z.weight_policy + '/baseline'
    z['coefficient_join_signed_difference'] = z.beta_dream-z.beta
    z['coefficient_join_difference'] = z.coefficient_join_signed_difference.abs()
    z['ICC_bin'] = pd.cut(z.ICC, ICC_EDGES, right=False, labels=ICC_LABELS)
    return z


def icc_summaries(z):
    summaries, strata = [], []
    for (spec, policy), g in z.groupby(['spec', 'weight_policy'], sort=True):
        available = g.status.eq('available')
        summaries.append(dict(spec=spec, weight_policy=policy, planned=len(g), available=int(available.sum()),
            failed=int((~available).sum()), nonconverged=int(g.converged.eq(False).sum()),
            singular_fraction=g.loc[available, 'singular'].mean(),
            zero_variance_fraction=g.loc[available, 'at_zero'].mean(),
            ICC_median=g.ICC.median(), ICC_q25=g.ICC.quantile(.25), ICC_q75=g.ICC.quantile(.75),
            ICC_q95=g.ICC.quantile(.95), within_cage_correlation_median=g.mean_within_cage_correlation.median()))
        for label in ICC_LABELS:
            a = g[g.ICC_bin.eq(label)]
            strata.append(dict(spec=spec, weight_policy=policy, ICC_bin=label, n=len(a),
                median_ICC=a.ICC.median(), median_SE_ratio=a.SE_ratio.median(),
                q25=a.SE_ratio.quantile(.25), q75=a.SE_ratio.quantile(.75)))
    return pd.DataFrame(summaries), pd.DataFrame(strata)


def correction_tables(human, mouse):
    h = human.rename(columns={'model': 'branch', 'contrast': 'model'}).copy()
    h.branch = h.branch.replace({'quality': 'extended'})
    h['collection'] = 'human_17064'
    sources = [h]
    for label, models in [('mouse_2844', ['baseline', 'cage_fixed']),
                          ('mouse_4266', ['baseline', 'cage_fixed', 'cage_correlation'])]:
        z = mouse[mouse.model.isin(models)].rename(columns={'model': 'branch', 'spec': 'model'}).copy()
        z.branch = z.branch.replace({'cage_correlation': 'cage_block'})
        z['collection'] = label
        z['human_gene_symbol'] = z.gene_symbol
        z['gene_symbol'] = z.mouse_gene.fillna('UNMAPPED:' + z.human_gene_symbol)
        sources.append(z)
    values, families, summaries, retention, transitions, continuous = [], [], [], [], [], []
    for z in sources:
        collection = z.collection.iloc[0]
        is_human = collection == 'human_17064'
        scopes = ({'original': ['branch'], 'by_contrast': ['branch', 'model'],
                   'by_context': ['branch', 'context'], 'by_context_contrast': ['branch', 'context', 'model']}
                  if is_human else {'original': [], 'by_method': ['branch'],
                   'by_specification': ['model'], 'by_specification_method': ['model', 'branch']})
        qcols = []
        for scope, fields in scopes.items():
            groups = list(z.groupby(fields, dropna=False, sort=True)) if fields else [('all', z)]
            for denominator in ['planned', 'available_only']:
                column = f'q__{scope}__{denominator}'
                qcols.append(column)
                z[column] = np.nan
                for family, g in groups:
                    p = g.p.to_numpy(float)
                    valid = np.isfinite(p)
                    q = reserved_bh(p) if denominator == 'planned' else np.full(len(p), np.nan)
                    if denominator == 'available_only' and valid.any():
                        q[valid] = reserved_bh(p[valid])
                    z.loc[g.index, column] = q
                    families.append(dict(collection=collection, scope=scope, family=str(family),
                        denominator=denominator, planned=len(g), available=int(valid.sum()),
                        unavailable=int((~valid).sum()), correction_n=len(g) if denominator=='planned' else int(valid.sum())))
        for (model, branch), g in z.groupby(['model', 'branch'], sort=True):
            original = g.q__original__planned.lt(.05)
            for scope in scopes:
                for denominator in ['planned', 'available_only']:
                    selected = g[f'q__{scope}__{denominator}'].lt(.05)
                    summaries.append(dict(collection=collection, model=model, branch=branch, scope=scope,
                        denominator=denominator, planned=len(g), available=int(g.p.notna().sum()),
                        significant=int(selected.sum()), original_significant=int(original.sum()),
                        added_vs_original=int((selected & ~original).sum()), lost_vs_original=int((original & ~selected).sum())))
        for model, g in z.groupby('model', sort=True):
            baseline = 'base' if is_human else 'baseline'
            keys = ['context', 'gene_symbol'] if is_human else ['gene_symbol']
            a = unique(g[g.branch.eq(baseline)], keys, 'BH baseline').set_index(keys)
            for branch, part in g[g.branch.ne(baseline)].groupby('branch', sort=True):
                b = unique(part, keys, 'BH comparator').set_index(keys).reindex(a.index)
                common = dict(collection=collection, model=model, comparison=branch)
                selected = a.q__original__planned.lt(.05)
                for col in qcols:
                    x, y = a[col], b[col]
                    both = x.notna() & y.notna()
                    retention.append(dict(**common, correction=col, original_base_discoveries=int(selected.sum()),
                        adjusted_significant=int(y[selected].lt(.05).sum()), adjusted_unavailable=int(y[selected].isna().sum())))
                    transitions.append(dict(**common, correction=col, planned=len(a), jointly_available=int(both.sum()),
                        base_only=int((x.lt(.05)&~y.lt(.05)&both).sum()), adjusted_only=int((~x.lt(.05)&y.lt(.05)&both).sum()),
                        both=int((x.lt(.05)&y.lt(.05)&both).sum()), neither=int((~x.lt(.05)&~y.lt(.05)&both).sum()),
                        unavailable=int((~both).sum())))
                t = a[[]].copy()
                for k, v in common.items(): t[k] = v
                for f in ['beta', 'SE', 'p']:
                    t[f+'_base'], t[f+'_extended'] = a[f], b[f]
                t['original_q_base'], t['original_q_extended'] = a.q__original__planned, b.q__original__planned
                with np.errstate(divide='ignore', invalid='ignore'):
                    t['log_abs_beta_ratio'] = np.log(abs(t.beta_extended/t.beta_base))
                    t['minus_log_SE_ratio'] = -np.log(t.SE_extended/t.SE_base)
                t['log_abs_t_ratio'] = t.log_abs_beta_ratio+t.minus_log_SE_ratio
                continuous.append(t.reset_index())
        values.append(z)
    return dict(values=pd.concat(values, ignore_index=True), families=pd.DataFrame(families),
        summary=pd.DataFrame(summaries), retention=pd.DataFrame(retention),
        transitions=pd.DataFrame(transitions), continuous=pd.concat(continuous, ignore_index=True))


def mouse_geometry(units):
    """Unweighted design geometry, not expression-model refitting."""
    def matrix(d, kind, cage=False):
        parts = [np.ones((len(d), 1)), d.genotype.eq('5X').to_numpy(float)[:, None]]
        parts.append((d.age-d.age.mean()).to_numpy(float)[:, None] if kind=='linear'
                     else pd.get_dummies(d.age, drop_first=True, dtype=float).to_numpy())
        parts.append(d.sex.eq('M').to_numpy(float)[:, None])
        if kind=='age_factor_slide':parts.append(pd.get_dummies(d.slide, drop_first=True, dtype=float).to_numpy())
        if cage:parts.append(pd.get_dummies(d.cage, drop_first=True, dtype=float).to_numpy())
        return np.column_stack(parts)
    mixed = units.groupby('cage').genotype.nunique()
    mixed_mask = units.cage.isin(mixed[mixed > 1].index).to_numpy()
    rows = []
    for region in ['whole', 'cortex']:
        for kind in ['linear', 'age_factor', 'age_factor_slide']:
            for branch, cage in [('baseline', False), ('cage_fixed', True), ('mixed_cages_geometry_only', True)]:
                x = matrix(units, kind, cage)
                if branch == 'mixed_cages_geometry_only': x = x[mixed_mask]
                r = x[:, 1]-np.delete(x, 1, axis=1) @ np.linalg.lstsq(np.delete(x, 1, axis=1), x[:, 1], rcond=None)[0]
                ss = float(r@r); rank = int(np.linalg.matrix_rank(x))
                rows.append(dict(spec=region+'_'+kind, branch=branch, n=len(x), columns=x.shape[1], rank=rank,
                    residual_df=len(x)-rank, genotype_residual_SS=ss,
                    unit_variance_beta_factor=1/ss if ss>1e-12 else np.nan))
    return pd.DataFrame(rows)


def precision_summaries(z):
    rows, errors, display = [], [], []
    for spec, g in z.groupby('spec', sort=True):
        rows.append(dict(spec=spec, genes=len(g), median_SE_ratio=g.SE_ratio.median(),
            median_unscaled_ratio=g.unscaled_ratio.median(), median_posterior_variance_ratio=g.posterior_variance_ratio.median(),
            median_shared_weights_SE_ratio=(g.SE_shared_weights/g.SE_base).median(), median_reweighting_SE_ratio=g.reweight_factor.median()))
        rmse = np.sqrt(np.mean(g.design_only_log_error**2))
        errors.append(dict(spec=spec, genes=len(g), design_only_ratio=g.design_only_ratio.iloc[0], observed_median=g.SE_ratio.median(),
            identity_max_log_error=np.abs(np.log(g.identity_ratio/g.SE_ratio)).max(),
            three_factor_max_log_error=np.abs(np.log(g.common_weight_identity_ratio/g.SE_ratio)).max(),
            design_only_median_log_error=g.design_only_log_error.median(), design_only_RMSE_log=rmse,
            design_only_q025_log_error=g.design_only_log_error.quantile(.025), design_only_q975_log_error=g.design_only_log_error.quantile(.975),
            weighted_information_RMSE_log=np.sqrt(np.mean(g.weighted_information_log_error**2))))
        display.append(dict(spec=spec, predicted=g.design_only_ratio.iloc[0], observed_median=g.SE_ratio.median(),
            observed_q25=g.SE_ratio.quantile(.25), observed_q75=g.SE_ratio.quantile(.75), log_RMSE=rmse))
    return pd.DataFrame(rows), pd.DataFrame(errors), pd.DataFrame(display)


def flat_json(value, path=''):
    if isinstance(value, dict):
        return [row for key, item in value.items() for row in flat_json(item, path+'/'+str(key))]
    if isinstance(value, list):
        return [dict(path=path, type='array', value=json.dumps(value, sort_keys=True))]
    return [dict(path=path, type=type(value).__name__, value=json.dumps(value))]


def clinical_change_summary(clinical):
    """All finite paired log components; selection is a separate count."""
    rows = []
    for context, g in clinical.groupby('context', sort=True):
        valid = np.isfinite(g[['log_abs_beta_ratio', 'minus_log_SE_ratio']]).all(axis=1)
        a = g.loc[valid]
        row = dict(context=context, n_testable=int(g.beta_base.notna().sum()),
                   n_decomposable=len(a), n_selected=int(g.original_q_base.lt(.05).sum()))
        for prefix, column in [('beta_log', 'log_abs_beta_ratio'), ('se_log', 'minus_log_SE_ratio')]:
            row[prefix+'_median'] = a[column].median()
            row[prefix+'_q25'] = a[column].quantile(.25)
            row[prefix+'_q75'] = a[column].quantile(.75)
        rows.append(row)
    return pd.DataFrame(rows)


def layout_membership(layout, five_model_results, three_model_results):
    """Membership comes from each experiment's results, not feasible geometry."""
    z = unique(layout.copy(), ['layout'], 'layout registry')
    for column, results in [('simulated_five_model', five_model_results),
                            ('simulated_three_model', three_model_results)]:
        names = set(results.dropna())
        require(names.issubset(set(z.layout)), 'Unregistered layouts in '+column)
        z[column] = z.layout.isin(names)
    z['simulated_here'] = z.simulated_five_model
    return z


def resource_catalog():
    roles = {
        'GSE118553':'Candidate-panel discovery input', 'SEA-AD':'Primary human and donor-conditioned analyses',
        'GSE233208':'Primary mouse analysis and observed design geometry', 'GSE48350':'Supplementary bulk endpoints',
        'GSE1297':'Supplementary cognition associations', 'GSE168137':'Supplementary external endpoints',
        'GSE157827':'Supplementary external endpoints', 'GSE92926':'Supplementary external endpoints',
        'GSE134151':'Supplementary external endpoints', 'GSE243292':'Astrocyte reference input',
        'GSE167490':'Supplementary astrocyte comparison', 'GSE222007':'Supplementary interactions and mixtures'}
    rows = []
    for accession, role in roles.items():
        sea = accession == 'SEA-AD'
        spec = ('DFC RNA-seq final-nuclei release 2026-06-22.' if sea else 'GEO accession '+accession+'.')
        rows.append(dict(resource=accession, analysis_role=role,
            relationship='primary' if accession in ('SEA-AD','GSE233208') else ('panel_selection' if accession=='GSE118553' else 'supplementary'),
            public_source='https://brain-map.org/consortia/sea-ad' if sea else 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+accession,
            input_specification=spec+' Input definitions and download code: public release v2.1.0.'))
    return pd.DataFrame(rows)


def simulation_jobs(protocol):
    """Expand the complete protocol array without changing job definitions."""
    jobs = protocol.get('jobs')
    require(isinstance(jobs,list) and jobs, 'Missing simulation job definitions')
    fields = ['id','phase','start','end','spec']
    require(all(isinstance(j,dict) and set(j)==set(fields) for j in jobs), 'Unexpected simulation job schema')
    require(len({j['id'] for j in jobs})==len(jobs), 'Duplicate simulation job IDs')
    return pd.DataFrame([dict(job_id=j['id'], phase=j['phase'], start=j['start'], end=j['end'],
        specification=json.dumps(j['spec'],sort_keys=True,ensure_ascii=False)) for j in jobs])


def complete(ex, root, supplements):
    """Rebind blocks to their declared statistical objects, not just a file count."""
    read = lambda p: pd.read_csv(p, low_memory=False)
    c, m, b = root/'core/results', root/'matched/results', root/'benchmark/final_results'
    def put(block, name, frame, definition):
        for row in ex.index:
            row['blocks'] = ','.join(x for x in row['blocks'].split(',') if x and x != block)
        return ex.add(f'S{int(block.split(".")[0][1:]):02}', name, frame=frame, blocks=block, description=definition)
    h, mouse = read(c/'human_results.csv.gz'), read(c/'mouse_results.csv.gz')
    put('S1.2','publication_resources',resource_catalog(),'Resource roles, public data entry points and input specifications.')
    t = correction_tables(h, mouse)
    for block, key in [('S3.7', 'families'), ('S3.8', 'retention'), ('S3.9', 'summary'),
                       ('S3.10', 'transitions'), ('S3.11', 'values'), ('S4.3', 'continuous')]:
        put(block, 'publication_'+key, t[key], 'Human and two separate mouse BH families; planned and available-only denominators.')
    summary = t['summary']; retention = t['retention']
    clinical = summary[summary.collection.eq('human_17064') & summary.model.eq('clinical_status')]
    retain = retention[retention.collection.eq('human_17064') & retention.model.eq('clinical_status')]
    put('S3.6', 'publication_family_display', pd.concat([clinical.assign(record_type='discoveries'), retain.assign(record_type='retention')], ignore_index=True), 'Clinical discovery and selected-subset retention views.')
    ch = t['continuous']; ch = ch[ch.collection.eq('human_17064') & ch.model.eq('clinical_status')]
    bins = []
    selected = ch[ch.original_q_base.lt(.05)]
    edges = [0, .01, .02, .03, .04, .05, .1, .2]
    for arm in ['base', 'extended']:
        q = selected['original_q_'+arm]
        for low, high in zip(edges[:-1], edges[1:]):
            bins.append(dict(arm=arm, lower=low, upper=high, n=int((q.ge(low)&q.lt(high)).sum())))
    put('S3.5', 'publication_q_bins', pd.DataFrame(bins), 'Selected original clinical discoveries, explicit half-open q intervals.')
    put('S4.2', 'publication_clinical_changes', clinical_change_summary(ch),
        'All finite paired clinical coefficient/SE log components by context; n_selected is a separate original-q<0.05 count.')
    put('S4.4', 'publication_selected_log_components', selected, 'All originally selected clinical associations, not all genes.')
    three = mouse[mouse.model.isin(['baseline', 'cage_fixed', 'cage_correlation'])]
    put('S12.1', 'publication_three_mouse_arms', three, 'Three-arm panel fits; seven-arm background results supplied separately.')
    from .reporting import pairs
    emp = pd.concat([pairs(h[h.gene_symbol.eq('EMP3')]).assign(resource='human'),
                     three[three.gene_symbol.eq('EMP3')].assign(resource='mouse')], ignore_index=True)
    put('S13.2', 'publication_EMP3', emp, 'EMP3: paired human contexts plus three mouse arms.')
    consistency = []
    for gene, g in three.groupby('mouse_gene'):
        valid = g.p.notna(); n = int(valid.sum())
        if n != 18: continue
        pos, neg = int(g.loc[valid, 'beta'].gt(0).sum()), int(g.loc[valid, 'beta'].lt(0).sum())
        consistency.append(dict(gene_symbol=gene, planned=len(g), testable=n, unavailable=len(g)-n,
            significant=int(g.q_three_model_joint.lt(.05).sum()), positive=pos, negative=neg,
            nonzero_same_direction_fraction=max(pos, neg)/(pos+neg) if pos+neg else np.nan))
    put('S13.3', 'publication_mouse_consistency', pd.DataFrame(consistency), 'Genes testable across all 18 three-arm/specification positions.')
    put('S13.4', 'publication_cortex_examples', three[three.spec.eq('cortex_age_factor') & three.gene_symbol.isin(['GFAP','GSDMD','IGFBP5'])], 'Cortex categorical-age specification, three genes by three arms.')
    z = native_precision(read(c/'mouse_SE_decomposition_check.csv.gz'), read(c/'mouse_precision.csv.gz'))
    ps, errors, display = precision_summaries(z)
    geo = mouse_geometry(read(root/'core/prepared/mouse_units.csv'))
    put('S14.1', 'publication_geometry', geo, 'Unweighted baseline, fixed and mixed-cage-only genotype information; not expression fits.')
    put('S14.2', 'publication_native_precision', z, 'Native fixed/base two-factor decomposition; separately named shared-weight factors.')
    put('S14.3', 'publication_precision_summary', ps, 'Medians of native factors and shared-weight contrast.')
    put('S14.5', 'publication_precision_plot_data', pd.concat([z.assign(record_type='gene'), geo.assign(record_type='geometry')], ignore_index=True), 'Gene decomposition and unweighted design geometry, distinguished by record_type.')
    put('S16.1', 'publication_design_errors', errors, 'Design-only prediction errors across genes; natural-log scale.')
    put('S16.2', 'publication_gene_design_errors', z[KEY+['observed_ratio','identity_ratio','common_weight_identity_ratio','design_only_ratio','design_only_log_error','weighted_information_log_error']], 'One row per gene and specification.')
    put('S16.3', 'publication_precision_definitions', pd.DataFrame(flat_json(dict(
        native_identity='SE_fixed/SE_base = unscaled_ratio * sqrt(posterior_variance_ratio)',
        shared_identity='SE_fixed/SE_base = design_factor * posterior_factor * reweight_factor',
        ICC_SE_ratio='dream_policy_SE / baseline_SE', ICC_edges=[.01,.1,.2],
        design_prediction='sqrt(residual_genotype_SS_base / residual_genotype_SS_fixed)', log_error='ln(observed_ratio/design_only_ratio)'))), 'Explicit numerical definitions; no archived-result comparison.')
    put('S16.4', 'publication_design_display', display, 'Design-only predicted ratio, observed median/IQR and log RMSE.')
    fits = mouse[mouse.model.eq('baseline')][KEY+['beta','SE']].rename(columns={'beta':'beta_base','SE':'SE_base'})
    mp = mouse.merge(fits, on=KEY, validate='many_to_one')
    mp['SE_ratio'] = mp.SE/mp.SE_base
    mp['absolute_change'] = (mp.beta-mp.beta_base).abs()
    mp['sign_flip'] = (mp.beta*mp.beta_base).lt(0)
    mp['q_spec'] = mp.groupby(['spec','model']).p.transform(reserved_bh)
    rows = []
    for (spec, model), g in mp.groupby(['spec','model']):
        rows.append(dict(spec=spec, model=model, planned=len(g), available=int(g.p.notna().sum()),
            significant_method_family=int(g.q_within_model.lt(.05).sum()), significant_spec_family=int(g.q_spec.lt(.05).sum()),
            median_SE_ratio=g.SE_ratio.median(), median_absolute_beta_change=g.absolute_change.median(), sign_flips=int(g.sign_flip.sum())))
    put('S14.4', 'publication_model_precision', pd.DataFrame(rows), 'Includes baseline and all five original arms; explicit method/specification BH.')
    fix = mp[mp.model.eq('cage_fixed')].copy()
    fix['base_significant'] = mouse[mouse.model.eq('baseline')].set_index(KEY).q_base_fixed_joint.reindex(pd.MultiIndex.from_frame(fix[KEY])).to_numpy() < .05
    fix['cage_significant'] = fix.q_base_fixed_joint.lt(.05)
    fs = fix.groupby('spec').agg(testable=('p','count'), base_significant=('base_significant','sum'), cage_significant=('cage_significant','sum'), median_SE_ratio=('SE_ratio','median'), median_absolute_change=('absolute_change','median')).reset_index()
    put('S13.1', 'publication_fixed_comparison', fs, 'Fixed versus baseline under the joint two-arm planned family.')
    cal = read(c/'calibration/calibration_summary.csv'); cal['SE_ratio'] = cal.RMS_SE/cal.empirical_SD
    agg = cal.groupby(['error_law','mean_case','method']).agg(
        rejection_median=('rejection','median'), rejection_min=('rejection','min'), rejection_max=('rejection','max'),
        coverage_median=('coverage','median'), coverage_min=('coverage','min'), coverage_max=('coverage','max'),
        SE_ratio_median=('SE_ratio','median'), power_or_rejection_median=('rejection','median'), joint_F_median=('joint_F_rejection','median')).reset_index()
    put('S20.4', 'publication_donor_calibration_summary', agg, 'Median and range across donor-conditioned designs, grouped by law/mean/method.')
    gen = read(b/'paired_summary.csv'); reg = gen[gen.status.eq('regular') & gen.mean_case.eq('quality_zero')]
    agg = reg.groupby(['error_law','method']).agg(regular_designs=('design','size'),
        rejection_min=('rejection','min'), rejection_median=('rejection','median'), rejection_max=('rejection','max'),
        coverage_min=('coverage','min'), coverage_median=('coverage','median'), coverage_max=('coverage','max')).reset_index()
    put('S20.8', 'publication_null_ranges', agg, 'Regular designs with quality_zero mean only; cancellation/power remain in full table.')
    failures = read(b/'cage_failures_and_warnings.csv.gz')
    put('S17.5', 'publication_numerical_failures', failures[failures.status.ne('not_estimable')], 'Numerical failures/warnings; structural unavailability retained in full failure export.')
    units = []
    geometry = read(b/'cage_geometry.csv')
    for layout in geometry.layout:
        units.append(read(b/'inputs'/(layout+'.csv')).assign(layout=layout))
    put('S19.6', 'publication_five_model_units', pd.concat(units, ignore_index=True), 'Exact saved inputs to the five-model benchmark, not the supplemental three-model simulation.')
    cage = read(b/'cage_summary.csv').merge(geometry, on='layout', validate='many_to_one')
    base = cage[cage.model.eq('baseline')][['scenario','mean_width']].rename(columns={'mean_width':'baseline_width'})
    cage = cage.merge(base, on='scenario', validate='many_to_one')
    cage['within_information_fraction'] = cage.SS_fixed/cage.SS_baseline
    cage['width_to_baseline'] = cage.mean_width/cage.baseline_width
    put('S19.3', 'publication_allocation_comparison', cage, 'Five-model results with matched baseline widths and layout geometry.')
    focal = read(supplements['astro']/'focal_results.csv')
    main = focal[focal.task.eq('astro_primary_all_donors')]
    require(len(main)==3, 'Missing primary three-gene astrocyte comparison')
    put('S8.3', 'publication_primary_astro', main, 'Primary AD/control astrocyte comparison; secondary family results remain in S8.1.')
    cal_geometry = read(c/'calibration/calibration_geometry.csv')
    cal_geometry['observed_df'] = cal_geometry['df']
    put('S7.4', 'publication_calibration_designs', cal_geometry, 'Actual donor-conditioned calibration matrices and residual df.')
    interactions = read(supplements['interactions']/'298_all_fixed_interaction_results.csv')
    put('S11.6', 'publication_filtered_interaction_discoveries',
        interactions[interactions.variant.eq('filtered_matched') & interactions.q_1422_within_variant.lt(.05)],
        'Descriptive display of filtered-branch discoveries; full family remains in S11.1.')
    layout = read(supplements['simulations']/'layout_registry.csv').rename(columns={
        'layout_id':'layout', 'mutant_fraction':'genotype_fraction', 'actual_mixed_fraction':'mixed_fraction',
        'genotype_SS_base':'SS_base','genotype_SS_fixed':'SS_fixed','requested_mutant_fraction':'requested_genotype_fraction'})
    layout['design_only_SE_ratio'] = np.sqrt(layout.SS_base/layout.SS_fixed.where(layout.SS_fixed>1e-12))
    put('S19.2', 'publication_layout_geometry', layout, 'Supplementary three-model layout geometry, including infeasible layouts.')
    three_model = read(supplements['simulations']/'model_summary.csv')
    put('S19.4', 'publication_layout_display', layout_membership(layout, cage.layout, three_model.layout_id),
        'Actual result-layout membership in five-model and supplementary three-model experiments; simulated_here is an alias of simulated_five_model, not feasibility.')
    for block, path in [('S17.7', b/'protocol.json'), ('S20.7', b/'protocol.json'),
                         ('S20.2', c/'calibration/calibration_plan.json'), ('S20.3', c/'calibration/calibration_status.json')]:
        if path.is_file():
            put(block, 'publication_'+block.replace('.','_')+'_settings', pd.DataFrame(flat_json(json.loads(path.read_text(encoding='utf-8')))), 'Current run configuration, flattened JSON paths; not manuscript equality evidence.')
    return z, geo
