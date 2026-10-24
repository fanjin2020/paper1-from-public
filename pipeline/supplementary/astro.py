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
from .common import save, digest, RDIR, need


FOCAL = ['IGFBP5', 'GSDMD', 'GFAP']


TARGET = 'GSE167490'


REFERENCE = 'GSE243292'


MIN_CELLS = 50


def geo(acc):
    return f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/'


def soft(path):
    entities, cur = [], None
    with gzip.open(path, 'rt', encoding='utf-8-sig') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('^'):
                kind, value = line[1:].split(' = ', 1)
                cur = dict(kind=kind, id=value, fields={})
                entities.append(cur)
            elif cur is not None and line.startswith('!') and ' = ' in line:
                key, value = line[1:].split(' = ', 1)
                cur['fields'].setdefault(key, []).append(value)
    return entities


def first(e, key, default=''):
    return e['fields'].get(key, [default])[0]


def chars(e):
    d = {}
    for value in e['fields'].get('Sample_characteristics_ch1', []):
        if ':' in value:
            key, val = value.split(':', 1)
            key, val = key.strip().lower(), val.strip()
            need(key not in d or d[key] == val, f'Conflicting metadata {e["id"]}/{key}')
            d[key] = val
    return d


def target_register(es):
    rows = []
    for e in es:
        if e['kind'] != 'SAMPLE':
            continue
        c = chars(e)
        m = re.fullmatch(r'Donor(\d+)(?:_replicate(\d+))?', first(e, 'Sample_title'))
        need(m is not None, f'Unrecognized author donor title: {e["id"]}')
        need(c.get('tissue') == 'Prefrontal cortex', 'Unexpected target brain region')
        disease = c.get('disease state', '')
        need(disease in ('Non-symptomatic', "Alzheimer's disease"), f'Unknown disease value: {disease!r}')
        sex = c.get('sex')
        need(sex in ('Male', 'Female'), 'Sex missing or unsupported')
        age = float(c['age'])
        need(18 <= age <= 120, 'Age out of range')
        links = [u.replace('ftp://', 'https://', 1) for k, vs in e['fields'].items()
                 if k.startswith('Sample_supplementary_file') for u in vs]
        row = dict(gsm=e['id'], donor='D' + m[1], replicate=m[2] or '1',
                   case=int(disease != 'Non-symptomatic'), age=age, sex=sex,
                   source_title=first(e, 'Sample_title'), disease_label=disease,
                   source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + e['id'])
        for name in ('barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz'):
            hits = [u for u in links if u.endswith(name)]
            need(len(hits) == 1, f'Expected one {name} for {e["id"]}')
            row[name] = hits[0]
        rows.append(row)
    frame = pd.DataFrame(rows)
    need(len(frame) == 17 and frame.gsm.is_unique and frame.donor.nunique() == 16, 'Target sample/donor inventory changed')
    for k in ('case', 'age', 'sex'):
        need((frame.groupby('donor')[k].nunique() == 1).all(), f'Technical replicates conflict: {k}')
    need(frame.drop_duplicates('donor').case.value_counts().to_dict() == {1: 10, 0: 6}, 'Target group counts changed')
    return frame


def decoded(a):
    return np.array([x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in a], dtype=object)


def h5_frame(group):
    result = {}
    for key, node in group.items():
        if key == '__categories':
            continue
        if isinstance(node, h5py.Group):
            need({'codes', 'categories'} <= set(node), f'Unsupported H5AD encoding: {key}')
            codes, cats = node['codes'][:], decoded(node['categories'][:])
            need(((codes >= -1) & (codes < len(cats))).all(), 'Invalid H5AD category code')
            result[key] = [cats[c] if c >= 0 else '' for c in codes]
        elif '__categories' in group and key in group['__categories']:
            codes, cats = node[:], decoded(group['__categories'][key][:])
            need(((codes >= -1) & (codes < len(cats))).all(), 'Invalid legacy category code')
            result[key] = [cats[c] if c >= 0 else '' for c in codes]
        else:
            a = node[:]
            result[key] = decoded(a) if a.dtype.kind in 'OSU' else a
    return pd.DataFrame(result)


def get_reference(dl):
    gz = dl.get(geo(REFERENCE) + 'suppl/GSE243292_ADsnRNAseq_GEO.h5ad.gz')
    p = dl.root / 'GSE243292_ADsnRNAseq_GEO.h5ad'
    rec = p.with_suffix('.decompressed.json')
    if p.exists():
        need(rec.exists() and readj(rec)['sha256'] == digest(p) and readj(rec)['compressed_sha256'] == digest(gz),
             'Reference decompression integrity mismatch')
    else:
        tmp = p.with_suffix('.part')
        with gzip.open(gz, 'rb') as src, tmp.open('wb') as dst:
            shutil.copyfileobj(src, dst, 8 * 1024**2)
        os.replace(tmp, p)
        writej(rec, dict(sha256=digest(p), compressed_sha256=digest(gz)))
    with h5py.File(p, 'r') as f:
        obs, var = h5_frame(f['obs']), h5_frame(f['var'])
        need({'sampleID', 'atscore', 'final_celltype', 'apoe', 'trem2'} <= set(obs), 'Reference obs schema changed')
        need('genename' in var, 'No reference gene symbols')
        node = f['X']
        enc = node.attrs.get('encoding-type', '')
        if isinstance(enc, bytes):
            enc = enc.decode()
        need(enc == 'csr_matrix', 'Only documented raw CSR X supported')
        shape = tuple(node.attrs['shape'])
        need(shape == (len(obs), len(var)), 'Reference dimensions disagree')
        x = sp.csr_matrix((node['data'][:], node['indices'][:], node['indptr'][:]), shape=shape)
    need(np.isfinite(x.data).all() and (x.data >= 0).all() and (x.data == np.floor(x.data)).all(), 'Reference X not counts')
    genes = var.genename.astype(str).to_numpy()
    unique, inv = np.unique(genes, return_inverse=True)
    x = x @ sp.csr_matrix((np.ones(len(inv), dtype=np.float32), (np.arange(len(inv)), inv)), shape=(len(inv), len(unique)))
    need(set(obs.atscore) == {'A-T-', 'A+T-', 'A+T+'}, 'Pathology coding changed')
    for field in ('atscore', 'apoe', 'trem2'):
        need((obs.groupby('sampleID')[field].nunique() == 1).all(), f'Reference donor inconsistency {field}')
    return x.tocsr(), unique, obs


def lognorm(x):
    lib = np.asarray(x.sum(axis=1)).ravel()
    need((lib > 0).all(), 'Zero count library')
    y = (sp.diags((10000 / lib).astype(np.float32)) @ x.astype(np.float32)).tocsr()
    y.data = np.log1p(y.data)
    return y


def classifier(x, genes, obs, forbidden, out, seed=271, allowed=None):
    """Reference-donor holdout, then locked full-reference annotation. No diagnosis used."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import confusion_matrix
    from threadpoolctl import threadpool_limits
    rng = np.random.default_rng(seed)
    labels = obs.final_celltype.astype(str).to_numpy()
    donors = obs.sampleID.astype(str).to_numpy()
    types = sorted(set(labels))
    need(len(types) == 7, 'Expected seven author reference types')
    astro = [c for c in types if c.lower() in ('ast', 'astro', 'astrocyte', 'astrocytes')]
    oligo = [c for c in types if c.lower() in ('oli', 'oligo', 'oligodendrocyte', 'oligodendrocytes')]
    need(len(astro) == len(oligo) == 1, 'Cannot identify explicit reference Astro/Oligo labels')
    available = np.array([g not in forbidden and (allowed is None or g in allowed)
                          and not g.startswith(('MT-', 'RPL', 'RPS')) for g in genes])
    need(available.sum() >= 500, 'Insufficient shared annotation gene universe')
    z = lognorm(x)
    test_donors = sorted(set(donors))[::3]
    train = ~np.isin(donors, test_donors)
    means = np.stack([np.asarray(z[train & (labels == c)].mean(axis=0)).ravel() for c in types])
    features = set()
    for j in range(len(types)):
        score = means[j] - np.max(np.delete(means, j, axis=0), axis=0)
        score[~available] = -np.inf
        features.update(np.argsort(-score, kind='stable')[:75])
    idx = np.array(sorted(features))
    need(len(idx) >= 100 and not set(genes[idx]) & forbidden, 'Annotation feature leakage')
    def sample_rows(mask):
        rows = []
        for d in sorted(set(donors[mask])):
            for c in types:
                ids = np.flatnonzero(mask & (donors == d) & (labels == c))
                rows.extend(rng.choice(ids, min(200, len(ids)), replace=False).tolist())
        return np.array(sorted(rows))
    def fit(rows):
        model = LogisticRegression(C=1.0, class_weight='balanced', max_iter=1500, solver='lbfgs', random_state=seed)
        with threadpool_limits(limits=2):
            model.fit(z[rows][:, idx], labels[rows])
        need((model.n_iter_ < 1500).all(), 'Annotation classifier did not converge')
        return model
    model = fit(sample_rows(train))
    held = np.flatnonzero(~train)
    prob = model.predict_proba(z[held][:, idx])
    pred, confidence, gap = classify_prob(model.classes_, prob)
    cm = confusion_matrix(labels[held], pred, labels=types + ['Uncertain'])
    save(out / '272_annotation_holdout_confusion.csv', pd.DataFrame(cm, index=types + ['Uncertain'], columns=types + ['Uncertain']).reset_index(names='true_label'))
    checks = []
    for c in types:
        tp = int(((pred == c) & (labels[held] == c)).sum())
        precision = tp / max(1, int((pred == c).sum()))
        recall = tp / max(1, int((labels[held] == c).sum()))
        checks.append(dict(cell_type=c, precision=precision, recall=recall,
                           n_true=int((labels[held] == c).sum()), n_predicted=int((pred == c).sum()),
                           gate_pass=precision >= .80 and recall >= .60))
    save(out / '272_annotation_holdout_metrics.csv', checks)
    save(out / '272_annotation_genes.csv', pd.DataFrame({'gene_symbol': genes[idx]}))
    writej(out / '272_annotation_scope.json', dict(reference=REFERENCE, held_out_donors=test_donors,
        method='Multinomial logistic reference label transfer; C=1, balanced classes',
        probability_threshold=.70, probability_gap=.20, excluded_all_fixed_endpoints=True,
        diagnosis_used_for_annotation=False, target_annotation_author_confirmed=False,
        target_accuracy_known=False, subtype_assignment='not_attempted',
        reference_holdout_does_not_prove_target_accuracy=True))
    gate = {r['cell_type']: bool(r['gate_pass']) for r in checks}
    need(gate[astro[0]], 'Astrocyte reference holdout gate failed; no target-label hypothesis testing')
    model = fit(sample_rows(np.ones(len(obs), dtype=bool)))
    return model, genes[idx], astro[0], oligo[0], gate


def classify_prob(classes, prob):
    order = np.argsort(prob, axis=1)
    confidence = prob[np.arange(len(prob)), order[:, -1]]
    gap = confidence - prob[np.arange(len(prob)), order[:, -2]]
    pred = np.asarray(classes, dtype=object)[order[:, -1]].copy()
    pred[(confidence < .70) | (gap < .20)] = 'Uncertain'
    return pred, confidence, gap


def read_library(dl, row):
    from threadpoolctl import threadpool_limits
    features = pd.read_csv(dl.get(row['features.tsv.gz']), sep='\t', header=None, dtype=str, keep_default_na=False)
    barcodes = pd.read_csv(dl.get(row['barcodes.tsv.gz']), sep='\t', header=None, dtype=str)[0].to_numpy()
    with gzip.open(dl.get(row['matrix.mtx.gz']), 'rb') as stream, threadpool_limits(limits=2):
        x = mmread(stream).tocsr().T.tocsr()
    x.sum_duplicates()
    need(x.shape == (len(barcodes), len(features)) and len(set(barcodes)) == len(barcodes), '10x dimensions/barcodes mismatch')
    need(features.shape[1] >= 2, 'Missing 10x gene symbol column')
    if features.shape[1] >= 3:
        keep = features[2].eq('Gene Expression').to_numpy()
        x, features = x[:, keep], features.loc[keep]
    need(np.isfinite(x.data).all() and (x.data >= 0).all() and (x.data == np.floor(x.data)).all(), 'Noninteger/negative target counts')
    symbols = features[1].to_numpy()
    genes, inv = np.unique(symbols, return_inverse=True)
    x = x @ sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(genes)))
    return x.tocsr(), genes, barcodes


def prepare_target(register, dl, model, annotation_genes, out):
    sums, census, qcs, annotations = {}, [], [], []
    common = None
    for row in register.to_dict('records'):
        print(f'272 processing {row["gsm"]} / {row["donor"]}', flush=True)
        x, genes, bc = read_library(dl, row)
        if common is None:
            common = genes
        need(np.array_equal(common, genes), 'Target libraries have different gene universes; do not silently intersect')
        loc = pd.Index(genes).get_indexer(annotation_genes)
        need((loc >= 0).all(), 'Missing locked annotation genes in target')
        lib = np.asarray(x.sum(axis=1)).ravel()
        detected = np.asarray((x > 0).sum(axis=1)).ravel()
        mt = np.asarray(x[:, np.char.startswith(genes.astype(str), 'MT-')].sum(axis=1)).ravel() / np.maximum(lib, 1)
        keep = (lib >= 500) & (lib <= 60000) & (detected >= 200) & (detected <= 7500) & (mt <= .05)
        labels = np.full(len(bc), 'QC_excluded', dtype=object)
        confidence, gap = np.full(len(bc), np.nan), np.full(len(bc), np.nan)
        if keep.any():
            prob = model.predict_proba(lognorm(x[keep])[:, loc])
            labels[keep], confidence[keep], gap[keep] = classify_prob(model.classes_, prob)
        save(out / (row['gsm'] + '_cell_annotations.csv.gz'), pd.DataFrame(dict(gsm=row['gsm'], donor=row['donor'],
             barcode=bc, label=labels, probability=confidence, margin=gap, total_umi=lib,
             detected_genes=detected, mitochondrial_fraction=mt, qc_pass=keep)))
        qcs.append(dict(gsm=row['gsm'], donor=row['donor'], cells=len(bc), qc_retained=int(keep.sum()),
                        uncertain=int((labels == 'Uncertain').sum()), doublet_removal='No dedicated doublet caller; residual doublets possible'))
        for cell in sorted(set(model.classes_) | {'Uncertain', 'AllRecovered'}):
            mask = keep if cell == 'AllRecovered' else labels == cell
            key = (row['donor'], cell)
            sums[key] = sums.get(key, np.zeros(len(genes))) + np.asarray(x[mask].sum(axis=0)).ravel()
            census.append(dict(gsm=row['gsm'], donor=row['donor'], cell_type=cell, n_cells=int(mask.sum())))
    counts = pd.DataFrame(census).groupby(['donor', 'cell_type'], as_index=False).n_cells.sum()
    save(out / '272_library_qc.csv', qcs)
    save(out / '272_library_cell_counts.csv', census)
    save(out / '272_donor_cell_counts.csv', counts)
    columns, meta = {}, []
    donors = register.drop_duplicates('donor').set_index('donor')
    for i, row in enumerate(counts.to_dict('records')):
        d, c = row['donor'], row['cell_type']
        uid = f'U{i:04d}'
        columns[uid] = sums[(d, c)]
        cov = donors.loc[d]
        meta.append(dict(unit_id=uid, donor=d, cell_type=c, n_cells=row['n_cells'],
                         case=int(cov['case']), age=float(cov.age), sex=cov.sex,
                         annotation='reference_transfer_not_author_label', biological_unit='donor'))
    x = pd.DataFrame(columns, index=common)
    m = pd.DataFrame(meta)
    x.to_csv(out / '272_pseudobulk_counts.csv.gz', index_label='gene_symbol')
    save(out / '272_analysis_units.csv', m)
    return x, m, counts


def fit(out, name, x, m, genes, formulas, rscript, mode='unpaired', coefficient='case'):
    task = out / name
    task.mkdir()
    need(m.unit_id.is_unique and set(m.unit_id) <= set(x.columns), 'Analysis unit mismatch')
    save(task / 'metadata.csv', m)
    save(task / 'endpoints.csv', pd.DataFrame({'gene_symbol': genes}))
    x.loc[:, m.unit_id].to_csv(task / 'expression.csv.gz', index_label='gene_symbol')
    (task / 'formulas.txt').write_text(mode + '\n' + ';'.join(formulas) + '\n' + coefficient + '\n', encoding='utf-8')
    with (task / 'engine.log').open('w', encoding='utf-8') as log:
        run = subprocess.run([rscript, '--vanilla', str(RDIR / 'astro.R'), str(task)], stdout=log, stderr=subprocess.STDOUT)
    need(run.returncode == 0, f'R engine failed: {task / "engine.log"}')
    z = pd.read_csv(task / 'endpoint_results.csv')
    z['task'] = name
    return z


def family(z, keys):
    z = z.copy()
    z['fixed_endpoint_BH'] = np.nan
    z['focal_BH'] = np.nan
    for _, idx in z.groupby(keys, dropna=False).groups.items():
        idx = list(idx)
        valid = z.loc[idx, 'P_Value'].notna()
        q = bh(z.loc[idx, 'P_Value'].fillna(1))
        z.loc[np.asarray(idx)[valid], 'fixed_endpoint_BH'] = q[valid]
        fi = [i for i in idx if z.loc[i, 'gene_symbol'] in FOCAL]
        need(len(fi) == 3, 'Focal family must keep all three slots')
        q = bh(z.loc[fi, 'P_Value'].fillna(1))
        valid = z.loc[fi, 'P_Value'].notna()
        z.loc[np.asarray(fi)[valid], 'focal_BH'] = q[valid]
    return z


def external_models(out, x, m, astro, fixed, rscript):
    z = m[(m.cell_type == astro) & (m.n_cells >= MIN_CELLS)].copy()
    need(z.donor.is_unique, 'Technical repeats not collapsed')
    z['age_c'] = z.age - z.age.mean()
    result = fit(out, 'astro_primary_all_donors', x, z, fixed, ['~age_c+factor(sex)+case'], rscript)
    # Author's published removal is sensitivity only; retain all eligible donors in primary.
    sub = z[~z.donor.isin(['D5', 'D9'])].copy()
    sub['age_c'] = sub.age - sub.age.mean()
    secondary = fit(out, 'astro_author_exclusion_sensitivity', x, sub, fixed, ['~age_c+factor(sex)+case'], rscript)
    result['scope'] = 'primary_external_followup_not_independence_certification'
    secondary['scope'] = 'secondary_author_donor_exclusion'
    joined = family(pd.concat([result, secondary], ignore_index=True), ['task', 'model'])
    save(out / '272_all_fixed_results.csv', joined)
    save(out / '272_focal_results.csv', joined[joined.gene_symbol.isin(FOCAL)])
    return joined


def ols_outcome(y, design):
    from scipy.stats import t
    y, a = np.asarray(y, float), np.asarray(design, float)
    need(np.isfinite(a).all() and np.isfinite(y).all(), 'Nonfinite OLS inputs')
    rank = np.linalg.matrix_rank(a)
    if rank < a.shape[1] or len(y) - rank < 4:
        return dict(testable=False, reason='rank_or_residual_df')
    b = np.linalg.lstsq(a, y, rcond=None)[0]
    df = len(y) - rank
    sigma = ((y - a @ b)**2).sum() / df
    r = np.linalg.qr(a, mode='reduced')[1]
    ri = np.linalg.solve(r, np.eye(r.shape[0]))
    se = np.sqrt(np.sum(ri * ri, axis=1) * sigma)
    if se[-1] <= 0:
        return dict(testable=False, reason='zero_residual_variance')
    crit = t.ppf(.975, df)
    return dict(testable=True, reason='tested', coefficient=float(b[-1]),
                CI_low=float(b[-1] - crit * se[-1]), CI_high=float(b[-1] + crit * se[-1]),
                P_Value=float(2 * t.sf(abs(b[-1] / se[-1]), df)), n=len(y), residual_df=int(df))


def composition_test(counts, donors, cell_types, out, prefix, adjustment):
    wide = counts.pivot(index='donor', columns='cell_type', values='n_cells').fillna(0)
    need('AllRecovered' in wide, 'Missing recovery denominator')
    m = donors.drop_duplicates('donor').set_index('donor').loc[wide.index]
    need(m.case.nunique() == 2, 'Composition comparison requires two groups')
    design = pd.DataFrame({'intercept': np.ones(len(m))}, index=m.index)
    if adjustment:
        design['age_c'] = m.age - m.age.mean()
        need(set(m.sex) <= {'Male', 'Female'}, 'Unknown sex')
        if m.sex.nunique() == 2:
            design['male'] = (m.sex == 'Male').astype(int)
    design['case'] = m.case.astype(float)
    results, fractions = [], []
    for c in cell_types:
        n = wide[c] if c in wide else pd.Series(0, index=wide.index)
        total = wide.AllRecovered
        need((n <= total).all(), 'Numerator exceeds denominator')
        valid = total >= MIN_CELLS
        y = np.log((n + .5) / (total - n + .5))
        if len(m.loc[valid]) < 8 or m.loc[valid, 'case'].value_counts().reindex([0, 1], fill_value=0).min() < 4:
            r = dict(testable=False, reason='insufficient_donors')
        else:
            r = ols_outcome(y[valid], design.loc[valid])
        results.append(dict(cell_type=c, scope=prefix, **r))
        for d in wide.index:
            fractions.append(dict(donor=d, case=int(m.loc[d, 'case']), cell_type=c, n_cells=int(n[d]),
                                  denominator=int(total[d]), fraction=float(n[d]/total[d]) if total[d] else None,
                                  logit_smoothed=float(y[d]), retained=bool(valid[d])))
    r = pd.DataFrame(results)
    if 'P_Value' not in r:
        r['P_Value'] = np.nan
    q = bh(r.P_Value.fillna(1))
    r['composition_BH'] = np.where(r.P_Value.notna(), q, np.nan)
    save(out / (prefix + '_proportion_results.csv'), r)
    save(out / (prefix + '_fractions.csv'), fractions)
    save(out / (prefix + '_design.csv'), design.reset_index())
    return pd.DataFrame(fractions)


def composition_and_state(out, x, m, census, obs, astro, oligo, gate, fixed, rscript):
    types = sorted(set(m.cell_type) - {'AllRecovered', 'Uncertain'})
    donor = m.drop_duplicates('donor')
    fractions = composition_test(census, donor, types + ['Uncertain'], out, '273_target_sorted_recovery', True)
    ref = obs[obs.atscore.isin(['A+T-', 'A+T+'])].copy()
    refcounts = ref.groupby(['sampleID', 'final_celltype']).size().rename('n_cells').reset_index()
    refcounts = refcounts.rename(columns={'sampleID': 'donor', 'final_celltype': 'cell_type'})
    totals = ref.groupby('sampleID').size().rename('n_cells').reset_index().rename(columns={'sampleID': 'donor'})
    totals['cell_type'] = 'AllRecovered'
    refcounts = pd.concat([refcounts, totals], ignore_index=True)
    refdonors = ref[['sampleID', 'atscore']].drop_duplicates().rename(columns={'sampleID': 'donor'})
    refdonors['case'] = (refdonors.atscore == 'A+T+').astype(int)
    composition_test(refcounts, refdonors, sorted(ref.final_celltype.unique()), out, '273_reference_author_recovery', False)
    mix = m[(m.cell_type == 'AllRecovered') & (m.n_cells >= MIN_CELLS)].copy()
    fr = fractions[fractions.cell_type == astro].set_index('donor')
    mix['astro_logit'] = mix.donor.map(fr.logit_smoothed)
    mix['age_c'] = mix.age - mix.age.mean()
    paired = fit(out, 'mixture_matched_composition_sensitivity', x, mix, fixed,
                 ['~age_c+factor(sex)+case', '~age_c+factor(sex)+astro_logit+case'], rscript)
    results = [paired]
    if gate.get(oligo, False):
        cm = m[m.cell_type.isin([astro, oligo]) & (m.n_cells >= MIN_CELLS)].copy()
        paired_donors = cm.groupby('donor').cell_type.nunique()
        cm = cm[cm.donor.isin(paired_donors[paired_donors == 2].index)].copy()
        cm['ct'] = np.where(cm.cell_type == astro, 'Astro', 'Oligo')
        cm['age_c'] = cm.age - cm.age.mean()
        results.append(fit(out, 'paired_astro_oligo_interaction', x, cm, fixed,
                           ['~ct*(age_c+factor(sex)+case)'], rscript,
                           mode='paired_celltype', coefficient='ctOligo:case'))
    else:
        writej(out / '273_interaction_not_testable.json', dict(reason='Oligo annotation holdout gate failed'))
        results.append(pd.DataFrame({'gene_symbol': list(fixed), 'task': 'paired_astro_oligo_interaction',
            'model': 1, 'logFC': np.nan, 'P_Value': np.nan, 'testable': False,
            'reason': 'Oligo annotation holdout gate failed'}))
    tab = family(pd.concat(results, ignore_index=True), ['task', 'model'])
    save(out / '273_all_state_context_results.csv', tab)
    save(out / '273_focal_state_context_results.csv', tab[tab.gene_symbol.isin(FOCAL)])
    wide = paired.pivot(index='gene_symbol', columns='model', values='logFC')
    if {1, 2} <= set(wide):
        wide.columns = ['base_coefficient', 'composition_adjusted_coefficient']
        wide['coefficient_difference_descriptive_only'] = wide.iloc[:, 1] - wide.iloc[:, 0]
        save(out / '273_matched_coefficient_comparison.csv', wide.reset_index())
    writej(out / '273_interpretation.json', dict(
        target_proportions='Conditional on LHX2-positive/NeuN-negative enrichment, QC and computational annotation; NOT tissue fractions',
        reference_proportions='Recovered nuclei, not stereological tissue abundance; reference donors reused, not new validation',
        mixture_adjustment='Conditional association, not mediation or causal removal of composition',
        interaction='Oligo minus Astro disease coefficient; donor block; not cross-species interaction',
        subtypes='Not inferred; no subtype-composition or subtype-causality claim',
        coefficient_difference='Descriptive only; no formal attenuation hypothesis test',
        independent_validation_certified=False))
