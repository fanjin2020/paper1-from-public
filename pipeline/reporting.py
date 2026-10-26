"""Export manuscript result families from one completed public-input experiment.

CSV exports are deliberately independent of Word/Excel layout. A block can use a
subset of a source table; source tables retain missing tests and full precision.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from scipy import stats

from .io import require, sha, dump
from .lifecycle import current_source, parent_run
from .supplementary.common import reserved_bh, discovery_direction_comparison
from .supplementary.driver import GROUPS
from . import publication_tables as publication


FOCAL = ('GFAP', 'GSDMD', 'IGFBP5')


def pairs(frame):
    keys = ['gene_symbol', 'spec', 'context', 'contrast']
    result = frame[frame.model.eq('base')].merge(frame[frame.model.eq('quality')],
        on=keys, suffixes=('_base', '_quality'), validate='one_to_one')
    result['delta'] = result.beta_quality - result.beta_base
    result['SE_ratio'] = result.SE_quality / result.SE_base
    result['sign_reversal'] = result.beta_base * result.beta_quality < 0
    result['base_significant'] = result.q_joint_base.lt(.05)
    result['quality_significant'] = result.q_joint_quality.lt(.05)
    return result


def scope_values(frame):
    """Keep one row per planned test and separate columns for each BH definition."""
    z = frame.copy()
    for scope, fields in {'all': [], 'contrast': ['contrast'], 'context': ['context'],
                          'context_contrast': ['spec']}.items():
        for denominator in ('planned', 'available'):
            def adjust(p):
                if denominator == 'planned':
                    return reserved_bh(p)
                p = np.asarray(p, float)
                valid = np.isfinite(p)
                q = np.full(len(p), np.nan)
                if valid.any():
                    q[valid] = reserved_bh(p[valid])
                return q
            z['q_' + scope + '_' + denominator] = z.groupby(['model'] + fields, dropna=False).p.transform(adjust)
    return z


def geometry(x):
    a = np.asarray(x, float)
    require(a.ndim == 2 and np.isfinite(a).all() and len(a), 'Invalid design matrix')
    rank = np.linalg.matrix_rank(a)
    variable = np.ptp(a, axis=0) > 0
    scaled = a.copy()
    scaled[:, variable] -= scaled[:, variable].mean(axis=0)
    norms = np.linalg.norm(scaled, axis=0)
    scaled[:, norms > 0] /= norms[norms > 0]
    row = dict(n=len(a), columns=a.shape[1], rank=int(rank), residual_df=len(a)-rank,
               scaled_condition_number=float(np.linalg.cond(scaled)))
    vifs = []
    for j in np.flatnonzero(variable):
        y = a[:, j]; other = np.delete(a, j, axis=1)
        residual = y - other @ np.linalg.lstsq(other, y, rcond=None)[0]
        total = np.sum((y-y.mean())**2); error = residual @ residual
        vifs.append(dict(predictor=x.columns[j], VIF=float(total/error) if error > 1e-20 else np.inf))
    return row, vifs


def donor_tables(metadata, repetitions=2000):
    keys = ['age', 'sex', 'pathology', 'PMI', 'RIN', 'pH', 'clinical_group']
    selected = metadata[metadata.n_cells.ge(20) & metadata.clinical_group.isin(['Dementia', 'No_dementia'])]
    selected = selected.dropna(subset=keys)
    require(selected.groupby('donor_id')[keys].nunique().le(1).all().all(), 'Inconsistent donor covariates')
    d = selected[['donor_id'] + keys].drop_duplicates().sort_values('donor_id').reset_index(drop=True)
    rng = np.random.default_rng(307)
    rows, associations = [], []
    for field in ['age', 'pathology', 'PMI', 'RIN', 'pH']:
        a = d.loc[d.clinical_group.eq('Dementia'), field].to_numpy(float)
        b = d.loc[d.clinical_group.eq('No_dementia'), field].to_numpy(float)
        require(min(len(a), len(b)) >= 2, 'Insufficient complete-case donors')
        denominator = np.sqrt((a.var(ddof=1) + b.var(ddof=1))/2)
        differences = np.array([rng.choice(a, len(a), replace=True).mean() -
                                rng.choice(b, len(b), replace=True).mean() for _ in range(repetitions)])
        row = dict(variable=field, Dementia_n=len(a), No_dementia_n=len(b),
            mean_difference=a.mean()-b.mean(),
            standardized_mean_difference=(a.mean()-b.mean())/denominator if denominator else np.nan,
            bootstrap_low=np.quantile(differences, .025), bootstrap_high=np.quantile(differences, .975))
        for label, values in [('Dementia', a), ('No_dementia', b)]:
            row.update({label+'_mean': values.mean(), label+'_SD': values.std(ddof=1),
                label+'_Q25': np.quantile(values, .25), label+'_median': np.median(values),
                label+'_Q75': np.quantile(values, .75)})
        rows.append(row)
        if field in ('PMI', 'RIN', 'pH'):
            u, p = stats.mannwhitneyu(a, b, alternative='two-sided', method='asymptotic')
            associations.append(dict(variable=field, contrast='clinical_status', statistic=u,
                effect=2*u/(len(a)*len(b))-1, effect_type='rank_biserial', p=p))
            rho, p = stats.spearmanr(d[field], d.pathology)
            associations.append(dict(variable=field, contrast='Braak', statistic=rho,
                effect=rho, effect_type='Spearman_rho', p=p))
    association = pd.DataFrame(associations)
    association['q_family6'] = reserved_bh(association.p)
    sex = pd.crosstab(d.sex, d.clinical_group).reset_index()
    return d, pd.DataFrame(rows), sex, association


def verify_sources(root):
    run_id = current_source(root / 'core', parent_run())
    for name in ('matched', 'benchmark', 'supplementary'):
        marker = json.loads((root / name / 'fresh_run.json').read_text(encoding='utf-8'))
        require(marker['run_id'] == run_id and marker['output'] == str((root/name).resolve()),
                'Reporting cannot mix experiments: ' + name)
    supplements = {}
    for name in GROUPS:
        receipt = json.loads((root / 'supplementary' / (name+'.complete.json')).read_text(encoding='utf-8'))
        directory = Path(receipt['directory']).resolve()
        require(directory.is_relative_to(root/'supplementary') and directory.is_dir(), 'Invalid stage directory')
        for relative, digest in receipt['files'].items():
            path = (root/'supplementary'/relative).resolve()
            require(path.is_relative_to(root/'supplementary') and path.is_file() and sha(path) == digest,
                    'Incomplete supplementary stage: ' + relative)
        supplements[name] = directory
    return run_id, supplements


def read(path):
    require(Path(path).is_file(), 'Missing computed result: ' + str(path))
    return pd.read_csv(path)


def row_count(path):
    opener = gzip.open if path.suffix == '.gz' else open
    previous_limit = csv.field_size_limit()
    # Serialized protocol job lists exceed csv's default 128-KiB field limit.
    csv.field_size_limit(max(previous_limit, 2**31-1))
    try:
        with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
            return max(0, sum(1 for _ in csv.reader(handle))-1)
    finally:
        csv.field_size_limit(previous_limit)


class Exports:
    def __init__(self, root):
        self.root = root
        self.index = []
        self.names = set()

    def add(self, section, name, frame=None, source=None, blocks='', description=''):
        require((frame is None) != (source is None), 'Specify exactly one export source')
        directory = self.root / section; directory.mkdir(parents=True, exist_ok=True)
        filename = name + ('.csv.gz' if source is not None and Path(source).suffix == '.gz' else '.csv')
        target = directory / filename
        require(target not in self.names, 'Duplicate export: ' + str(target))
        self.names.add(target)
        if source is not None:
            require(Path(source).is_file(), 'Missing computed result: ' + str(source))
            shutil.copyfile(source, target)
        else:
            frame.to_csv(target, index=False)
        self.index.append(dict(section=section, blocks=blocks, file=target.relative_to(self.root).as_posix(),
            rows=row_count(target), bytes=target.stat().st_size, sha256=sha(target), description=description))
        return target

    def finish(self):
        index = pd.DataFrame(self.index)
        index.to_csv(self.root/'source_index.csv', index=False)
        rows = []
        for record in self.index:
            for block in record['blocks'].split(','):
                if block:
                    rows.append(dict(block=block, file=record['file'], description=record['description']))
        require(len({r['block'] for r in rows}) == len(publication.EXPECTED_BLOCKS), 'Incomplete workbook-block source map')
        require(len(rows) == len({r['block'] for r in rows}), 'A block must have one explicit combined source')
        expected=publication.EXPECTED_BLOCKS
        require({r['block'] for r in rows} == expected, 'Unexpected workbook block identifiers')
        pd.DataFrame(rows).to_csv(self.root/'block_sources.csv', index=False)
        from .workbook import generate as write_workbook
        write_workbook(self.root)
        return index


def generate(root, output, supplements):
    c = root/'core/results'; m = root/'matched/results'; b = root/'benchmark/final_results'
    p = root/'core/prepared'
    ex = Exports(output)
    h = read(c/'human_results.csv.gz'); hp = pairs(h)
    clinical = hp[hp.contrast.eq('clinical_status')]
    selected = clinical[clinical.base_significant]
    mouse = read(c/'mouse_results.csv.gz')
    changes = read(c/'human_paired_changes.csv.gz')

    ex.add('S01', 'cell_class_readiness', source=supplements['interactions']/'cell_class_readiness.csv', blocks='S1.1')
    accessions = ['GSE118553','SEA-AD','GSE233208','GSE48350','GSE1297','GSE168137','GSE157827',
                  'GSE92926','GSE134151','GSE243292','GSE167490','GSE222007']
    ex.add('S01', 'resources', frame=pd.DataFrame({'resource': accessions}), blocks='S1.2',
           description='Download URLs, sizes and SHA256 are in raw download receipts.')
    ex.add('S01', 'human_metadata', source=p/'human/metadata.csv', blocks='S1.3')
    ex.add('S02', 'human_paired_coefficients', frame=hp, blocks='S2.1')
    ex.add('S02', 'human_gene_membership', source=c/'human_gene_membership.csv', blocks='S2.2')
    checks = h.assign(CI_midpoint_error=(h.low+h.high)/2-h.beta,
                      CI_halfwidth_error=(h.high-h.low)/2-stats.t.ppf(.975,h.df)*h.SE)
    ex.add('S02', 'human_model_checks', frame=checks, blocks='S2.3',
           description='Within-run coefficient, SE, df and CI identities; no archived-fit comparison.')
    magnitude = []
    for contrast, z in hp.groupby('contrast'):
        for label, a in [('all_testable', z.dropna(subset=['beta_base','beta_quality'])),
                         ('sign_reversal', z[z.sign_reversal]), ('base_significant', z[z.base_significant])]:
            magnitude.append(dict(contrast=contrast, subset=label, n=len(a),
                median_absolute_base=a.beta_base.abs().median(), median_absolute_quality=a.beta_quality.abs().median(),
                median_absolute_change=a.delta.abs().median(), median_SE_ratio=a.SE_ratio.median()))
    ex.add('S03', 'coefficient_magnitudes', frame=pd.DataFrame(magnitude), blocks='S3.1')
    ex.add('S03', 'selected_clinical_associations', frame=selected, blocks='S3.2')
    denominator = clinical.groupby('context').agg(planned=('gene_symbol','size'),
        base_available=('p_base','count'), quality_available=('p_quality','count'),
        base_significant=('base_significant','sum'), quality_significant=('quality_significant','sum')).reset_index()
    ex.add('S03', 'clinical_contexts', frame=denominator, blocks='S3.3')
    ex.add('S03', 'clinical_pairs', frame=clinical, blocks='S3.4')
    bins = []
    for label, z in [('all',clinical),('base_significant',selected)]:
        for arm in ('base','quality'):
            q = z['q_joint_'+arm]
            cat = pd.cut(q, [0,.01,.05,.1,.2,.5,1], include_lowest=True, right=False).astype('string')
            cat[q.eq(1)] = '[0.5, 1.0)'; cat[q.isna()] = 'unavailable'
            for interval, n in cat.value_counts(dropna=False).items():
                bins.append(dict(subset=label,model=arm,q_interval=interval,n=n))
    ex.add('S03', 'q_bins', frame=pd.DataFrame(bins), blocks='S3.5')
    ex.add('S03', 'correction_scope_summary', source=c/'human_correction_scope_summary.csv', blocks='S3.6,S3.8,S3.9,S3.10',
           description='Selected discoveries and retention for each declared family; denominators are explicit.')
    ex.add('S03', 'correction_families', source=c/'correction_scope.csv', blocks='S3.7')
    ex.add('S03', 'correction_values', frame=scope_values(h), blocks='S3.11')
    ex.add('S04', 'paired_changes', frame=changes, blocks='S4.1')
    z = clinical.merge(changes, on=['gene_symbol','spec'], validate='one_to_one', suffixes=('','_paired'))
    z['change_model_significant'] = z.q_joint.lt(.05); z['change_HC3_significant'] = z.q_HC3_joint.lt(.05)
    summary = z.groupby('context').agg(planned=('gene_symbol','size'),
        changed_model=('change_model_significant','sum'),changed_HC3=('change_HC3_significant','sum'),
        sign_reversals=('sign_reversal','sum')).reset_index()
    ex.add('S04', 'clinical_continuous_changes', frame=summary, blocks='S4.2')
    ex.add('S04', 'continuous_transitions', frame=hp, blocks='S4.3', description='Human contrasts; mouse transitions in S13.')
    ex.add('S04', 'human_precision', source=c/'human_precision.csv.gz', blocks='S4.4')

    ex.add('S05', 'age_sensitivity', source=supplements['age']/'183_fixed_endpoint_summary.csv', blocks='S5.1')
    donors, balance, sexes, association = donor_tables(read(p/'human/metadata.csv'))
    for name, data, block in [('quality_balance',balance,'S5.2'),('sex_counts',sexes,'S5.3'),
                              ('quality_associations',association,'S5.4'),('donors',donors,'S5.5')]:
        ex.add('S05', name, frame=data, blocks=block)
    ex.add('S06','quality_model_variants',source=supplements['quality']/'quality_model_variants.csv',blocks='S6.1')
    ex.add('S06','quality_variant_genes',source=supplements['quality']/'quality_variants.csv')
    for i, name in enumerate(['human_sensitivity_discoveries','human_paired_family_and_selected_subset',
                              'human_donor_sex_counts','human_input_manifest'],2):
        ex.add('S06',name,source=m/(name+'.csv'),blocks='S6.'+str(i))
    ex.add('S06','human_context_sex_counts',source=m/'human_context_sex_counts.csv',blocks='S6.6',
           description='Human sex counts by clinical group and analysis context.')
    geo, vif = [], []
    for i, record in enumerate(read(c/'human_designs.csv').itertuples(index=False),1):
        for arm in ('base','quality'):
            path = p/'designs'/f'{i}_{arm}.csv'
            if path.is_file():
                x = read(path).drop(columns='unit_id'); row, v = geometry(x)
                geo.append(dict(spec=record.spec,model=arm,status='available',**row))
                vif.extend(dict(spec=record.spec,model=arm,**a) for a in v)
            else:
                require(record.status != 'available','Missing estimable design '+str(path))
                geo.append(dict(spec=record.spec,model=arm,status='not_estimable',n=record.n))
    ex.add('S07','design_geometry',frame=pd.DataFrame(geo),blocks='S7.1,S7.4',description='Geometry of this run actual matrices.')
    ex.add('S07','predictor_VIF',frame=pd.DataFrame(vif),blocks='S7.2')
    ex.add('S07','calibration_geometry',source=c/'calibration/calibration_geometry.csv',blocks='S7.3')
    ex.add('S08','astrocyte_candidates',source=supplements['astro']/'focal_results.csv',blocks='S8.1,S8.3')
    ex.add('S08','SEA_astrocyte_candidates',frame=h[h.context.eq('Astrocyte') & h.gene_symbol.isin(FOCAL)],blocks='S8.2')
    ex.add('S09','bulk_endpoints',source=supplements['bulk']/'188_all_fixed_endpoint_results.csv',blocks='S9.1')
    external = read(supplements['external']/'all_fixed_endpoint_results.csv')
    external['q_primary_focus_36'] = external.q_across_primary_directions
    external = discovery_direction_comparison(external)
    ex.add('S09','external_endpoints',frame=external,blocks='S9.2')
    dictionary = pd.DataFrame([
        ('P_Value','Two-sided within-model gene test before multiplicity adjustment'),
        ('logFC','Fitted coefficient on the declared cohort scale; consult formula and coefficient to distinguish main effects from interactions'),
        ('primary_task','Prespecified primary comparison, not independence certification'),
        ('primary_focus_36_member','Membership in the reserved 36-slot focal family'),
        ('primary_family_size','All planned slots, including untestable slots'),
        ('q_primary_focus_36','BH across the reserved primary 36-slot family; alias of q_across_primary_directions'),
        ('same_as_discovery_direction','Main-effect sign comparison, not replication evidence; missing for interactions, unavailable estimates or unknown discovery directions'),
        ('direction_comparison_status','comparable, not_applicable_interaction, unavailable_estimate or unknown_discovery_direction; does not change test availability')],columns=['field','definition'])
    ex.add('S09','field_dictionary',frame=dictionary,blocks='S9.3')
    ex.add('S09','primary_focal_results',frame=external[external.primary_focus_36_member.eq(True)],blocks='S9.4')
    ex.add('S10','cognition',source=supplements['cognition']/'231_all_fixed_gene_results.csv',blocks='S10.1')
    ex.add('S10','cognition_signatures',source=supplements['cognition']/'232_all_fixed_signature_results.csv',blocks='S10.2',
           description='CAMERA cognition gene-set tests.')
    interactions = read(supplements['interactions']/'298_all_fixed_interaction_results.csv')
    ex.add('S11','interactions',frame=interactions,blocks='S11.1,S11.6',description='Full family, including nonsignificant and missing tests.')
    ex.add('S11','doublet_seed_stability',source=supplements['interactions']/'doublet_seed_summary.csv',blocks='S11.2')
    ex.add('S11','doublet_libraries',source=supplements['interactions']/'doublet_library_summary.csv')
    mixture = read(supplements['mixtures']/'323_matched_cell_and_synthetic_mixture_results.csv')
    composition = read(supplements['mixtures']/'323_recovered_composition_associations.csv')
    ex.add('S11','cell_mixtures',frame=mixture,blocks='S11.3')
    ex.add('S11','composition_associations',frame=composition,blocks='S11.4')
    ex.add('S11','mixture_plot_data',frame=pd.concat([mixture.assign(record_type='expression'),
        composition.assign(record_type='composition')],ignore_index=True),blocks='S11.5')
    ex.add('S12','mouse_panel_models',frame=mouse,blocks='S12.1,S12.5')
    for i,name in [(2,'mouse_dream_diagnostics'),(3,'mouse_dream_retry_attempts'),(4,'mouse_model_status')]:
        ex.add('S12',name,source=c/(name+'.csv'),blocks='S12.'+str(i))
    keys = ['gene_symbol','mouse_gene','spec']
    baseline = mouse[mouse.model.eq('baseline')].drop(columns='model')
    mp = mouse[~mouse.model.eq('baseline')].merge(baseline,on=keys,suffixes=('','_base'),validate='many_to_one')
    mp['SE_ratio'] = mp.SE/mp.SE_base; mp['beta_difference'] = mp.beta-mp.beta_base
    mp['sign_reversal'] = mp.beta*mp.beta_base < 0
    ex.add('S12','mouse_paired_SE_ratios',frame=mp,blocks='S12.6')
    for i,name in enumerate(['mouse_method_summary','mouse_matched_method_summary','mouse_augmented_designs',
                            'weighted_centering_identity_summary','mouse_sex_by_genotype'],7):
        ex.add('S12',name,source=m/(name+'.csv'),blocks='S12.'+str(i))
    ms = mp.groupby(['spec','model']).agg(planned=('gene_symbol','size'),available=('p','count'),
            median_SE_ratio=('SE_ratio','median'),sign_reversals=('sign_reversal','sum')).reset_index()
    ex.add('S13','mouse_paired_summary',frame=ms,blocks='S13.1')
    ex.add('S13','EMP3_context',frame=mouse[mouse.gene_symbol.eq('EMP3')],blocks='S13.2')
    ex.add('S13','mouse_continuous_transitions',frame=mp,blocks='S13.3')
    ex.add('S13','mouse_cortex_examples',frame=mouse[mouse.spec.str.startswith('cortex') & mouse.gene_symbol.isin(FOCAL)],blocks='S13.4')
    precision = read(c/'mouse_precision.csv.gz')
    ps = precision.groupby('spec')[['SE_ratio','design_factor','posterior_factor','reweight_factor',
        'unweighted_design_prediction']].agg(['median','min','max']).reset_index()
    ps.columns = ['_'.join(x).rstrip('_') for x in ps.columns]
    ex.add('S14','mouse_geometry',source=c/'mouse_designs.csv',blocks='S14.1,S14.5')
    ex.add('S14','SE_components',frame=precision,blocks='S14.2,S14.5')
    ex.add('S14','SE_summary',frame=ps,blocks='S14.3')
    ex.add('S14','mouse_precision_by_model',frame=ms,blocks='S14.4')
    variance = read(c/'mouse_variance.csv.gz')
    iv = publication.icc_precision(variance, read(m/'mouse_seven_arm_background_fits.csv.gz'))
    ex.add('S15','ICC_precision',frame=iv,blocks='S15.1')
    ex.add('S15','variance_components',frame=variance,blocks='S15.2')
    icc_summary, icc_strata = publication.icc_summaries(iv)
    for name,tab,block in [('ICC_summary',icc_summary,'S15.3'),('ICC_strata',icc_strata,'S15.4')]:
        ex.add('S15',name,frame=tab,blocks=block)
    ex.add('S16','precision_summary',frame=ps,blocks='S16.1,S16.4')
    ex.add('S16','precision_gene_diagnostics',source=c/'mouse_SE_decomposition_check.csv.gz',blocks='S16.2')
    ex.add('S16','precision_design_plan',source=c/'mouse_designs.csv',blocks='S16.3')
    for i,name in [(1,'cage_summary'),(2,'cage_scenarios'),(3,'cage_paired_method_differences'),(4,'mundlak_fixed_identity')]:
        ex.add('S17',name,source=b/(name+'.csv'),blocks='S17.'+str(i))
    ex.add('S17','fit_failures_and_warnings',source=b/'cage_failures_and_warnings.csv.gz',blocks='S17.5')
    cage = read(b/'cage_summary.csv')
    ex.add('S17','model_availability',frame=cage.groupby('model')[['planned','available','not_estimable','numerical_failed']].sum().reset_index(),blocks='S17.6')
    protocol = json.loads((b/'protocol.json').read_text(encoding='utf-8'))
    settings = pd.DataFrame([dict(setting=k,value=json.dumps(v,sort_keys=True)) for k,v in protocol.items()])
    ex.add('S17','run_settings',frame=settings,blocks='S17.7')
    for i,name in enumerate(['322_design_simulation','model_summary','paired_model_comparisons'],1):
        ex.add('S18',name,source=supplements['simulations']/(name+'.csv'),blocks='S18.'+str(i))
    ex.add('S19','layout_registry',source=supplements['simulations']/'layout_registry.csv',blocks='S19.1,S19.2,S19.4')
    ex.add('S19','five_model_allocation',frame=cage,blocks='S19.3')
    ex.add('S19','five_model_geometry',source=b/'cage_geometry.csv',blocks='S19.5')
    ex.add('S19','layout_units',source=supplements['simulations']/'layout_units.csv',blocks='S19.6')
    calibration = read(c/'calibration/calibration_summary.csv')
    ex.add('S20','donor_conditioned_calibration',frame=calibration,blocks='S20.1,S20.4',description='All contexts, not just display averages.')
    for i,name in [(2,'plan'),(3,'status')]:
        path = c/'calibration'/(name+'.json')
        if not path.is_file():
            path = c/'calibration'/('calibration_'+name+'.json')
        obj = json.loads(path.read_text(encoding='utf-8'))
        ex.add('S20','calibration_'+name,frame=pd.DataFrame([dict(setting=k,value=json.dumps(v,sort_keys=True)) for k,v in obj.items()]),blocks='S20.'+str(i))
    general = read(b/'paired_summary.csv')
    ex.add('S20','general_calibration',frame=general,blocks='S20.5')
    ex.add('S20','general_geometry',source=b/'paired_geometry.csv',blocks='S20.6')
    ex.add('S20','general_settings',frame=settings,blocks='S20.7')
    reg = general[general.status.eq('regular')]
    ranges = reg.groupby(['error_law','mean_case','method'])[['rejection','coverage']].agg(['min','max']).reset_index()
    ranges.columns = ['_'.join(x).rstrip('_') for x in ranges.columns]
    ex.add('S20','general_ranges',frame=ranges,blocks='S20.8')
    for i,name in enumerate(['count_scenarios','count_pipeline_summary','count_pipeline_paired_summary','count_pipeline_diagnostic_inventory'],9):
        ex.add('S20',name,source=m/(name+'.csv'),blocks='S20.'+str(i))
    for i,name in enumerate(['count_pipeline_per_dataset','count_pipeline_paired_per_dataset',
                             'count_pipeline_unavailable_genes','count_generator_parameters'],13):
        ex.add('S20',name,source=m/(name+'.csv'),blocks='S20.'+str(i))
    ex.add('S20','simulation_jobs',frame=publication.simulation_jobs(protocol),blocks='S20.17',
           description='Complete fixed-design simulation job definitions; one job per row.')
    publication.complete(ex, root, supplements)
    return ex.finish()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args(argv)
    root = args.run.resolve()
    run_id, supplements = verify_sources(root)
    output = root/'reporting'; output.mkdir(exist_ok=True)
    dump(output/'status.json',dict(status='RUNNING',run_id=run_id))
    try:
        index = generate(root,output,supplements)
        from .supplementary_figures import generate as draw
        figures = draw(output)
        dump(output/'status.json',dict(status='RESULT_EXPORTS_COMPLETE',run_id=run_id,
             files=len(index), sections=20, workbook_blocks_mapped=len(publication.EXPECTED_BLOCKS), supplementary_figures=len(figures),
             workbook='Supplementary_Tables.xlsx', workbook_sha256=sha(output/'Supplementary_Tables.xlsx'),
             reporting_version=publication.VERSION,
             scope='Numeric source tables and 112-block workbook; manuscript agreement requires separate verification'))
    except BaseException as exc:
        dump(output/'status.json',dict(status='FAILED',run_id=run_id,reason=str(exc)))
        raise
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
