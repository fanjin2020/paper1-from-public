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
from .common import Blocked, utc, need
TARGET = 'GSE222007'
REFERENCE = 'GSE243292'
FOCAL = ['GFAP', 'GSDMD', 'IGFBP5']
MARKERS = {'astro': ['ALDH1L1', 'AQP4', 'SLC1A2', 'SLC1A3', 'GLUL', 'SOX9'], 'vascular': ['CLDN5', 'PECAM1', 'VWF', 'RGS5', 'PDGFRB', 'ACTA2', 'COL1A1'], 'oligo': ['PLP1', 'MOG', 'MBP']}
POLICY = dict(min_cells=50, min_group=8, min_df=8, min_age_range_retention=0.7, max_abs_age_smd=0.5, min_age_bin_overlap=0.5, min_umi=500, max_umi=60000, min_features=200, max_features=7500, max_mito=0.05, min_marker_genes=3, marker_log_difference=0.25, max_contaminant_fraction=0.2)
HGNC_URLS = ('https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt', 'https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt')

def geo(acc):
    need(re.fullmatch('GSE\\d+', acc), 'Bad GEO accession')
    return f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/'

def soft(path):
    entities, cur = ([], None)
    with gzip.open(path, 'rt', encoding='utf-8-sig') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('^'):
                kind, value = line[1:].split(' = ', 1)
                cur = dict(kind=kind, id=value, fields={})
                entities.append(cur)
            elif cur is not None and line.startswith('!') and (' = ' in line):
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
            key, val = (key.strip().lower(), val.strip())
            need(key not in d or d[key] == val, f'Conflicting metadata {e['id']}/{key}')
            d[key] = val
    return d

def decoded(a):
    return np.array([x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in a], dtype=object)

def h5_frame(group):
    result = {}
    for key, node in group.items():
        if key == '__categories':
            continue
        if isinstance(node, h5py.Group):
            need({'codes', 'categories'} <= set(node), f'Unsupported H5AD encoding: {key}')
            codes, cats = (node['codes'][:], decoded(node['categories'][:]))
            need(((codes >= -1) & (codes < len(cats))).all(), 'Invalid H5AD category code')
            result[key] = [cats[c] if c >= 0 else '' for c in codes]
        elif '__categories' in group and key in group['__categories']:
            codes, cats = (node[:], decoded(group['__categories'][key][:]))
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
        need(rec.exists() and readj(rec)['sha256'] == digest(p) and (readj(rec)['compressed_sha256'] == digest(gz)), 'Reference decompression integrity mismatch')
    else:
        tmp = p.with_suffix('.part')
        with gzip.open(gz, 'rb') as src, tmp.open('wb') as dst:
            shutil.copyfileobj(src, dst, 8 * 1024 ** 2)
        os.replace(tmp, p)
        writej(rec, dict(sha256=digest(p), compressed_sha256=digest(gz)))
    with h5py.File(p, 'r') as f:
        obs, var = (h5_frame(f['obs']), h5_frame(f['var']))
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
    need((x.data < 2 ** 53).all(), 'Reference counts exceed exact conversion range')
    genes = var.genename.astype(str).to_numpy()
    unique, inv = np.unique(genes, return_inverse=True)
    x = x.astype(np.int64) @ sp.csr_matrix((np.ones(len(inv), dtype=np.int64), (np.arange(len(inv)), inv)), shape=(len(inv), len(unique)))
    need(set(obs.atscore) == {'A-T-', 'A+T-', 'A+T+'}, 'Pathology coding changed')
    for field in ('atscore', 'apoe', 'trem2'):
        need((obs.groupby('sampleID')[field].nunique() == 1).all(), f'Reference donor inconsistency {field}')
    return (x.tocsr(), unique, obs)

def lognorm(x):
    lib = np.asarray(x.sum(axis=1)).ravel()
    need((lib > 0).all(), 'Zero count library')
    y = (sp.diags((10000 / lib).astype(np.float32)) @ x.astype(np.float32)).tocsr()
    y.data = np.log1p(y.data)
    return y

def classifier(x, genes, obs, forbidden, out, seed=271, allowed=None):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import confusion_matrix
    from threadpoolctl import threadpool_limits
    """Reference-donor holdout, then locked full-reference annotation. No diagnosis used."""
    rng = np.random.default_rng(seed)
    labels = obs.final_celltype.astype(str).to_numpy()
    donors = obs.sampleID.astype(str).to_numpy()
    need(obs.sampleID.notna().all() and all((d.strip() for d in donors)), 'Missing reference donor identity')
    types = sorted(set(labels))
    need(len(types) == 7, 'Expected seven author reference types')
    astro = [c for c in types if c.lower() in ('ast', 'astro', 'astrocyte', 'astrocytes')]
    oligo = [c for c in types if c.lower() in ('oli', 'oligo', 'oligodendrocyte', 'oligodendrocytes')]
    need(len(astro) == len(oligo) == 1, 'Cannot identify explicit reference Astro/Oligo labels')
    available = np.array([g not in forbidden and (allowed is None or g in allowed) and (not g.startswith(('MT-', 'RPL', 'RPS'))) for g in genes])
    need(available.sum() >= 500, 'Insufficient shared annotation gene universe')
    z = lognorm(x)
    test_donors = sorted(set(donors))[::3]
    train = ~np.isin(donors, test_donors)
    need(all(((train & (labels == c)).any() and (~train & (labels == c)).any() for c in types)), 'Reference donor split lacks at least one class in training or holdout')
    means = np.stack([np.asarray(z[train & (labels == c)].mean(axis=0)).ravel() for c in types])
    features = set()
    for j in range(len(types)):
        score = means[j] - np.max(np.delete(means, j, axis=0), axis=0)
        score[~available] = -np.inf
        features.update(np.argsort(-score, kind='stable')[:75])
    idx = np.array(sorted(features))
    need(len(idx) >= 100 and (not set(genes[idx]) & forbidden), 'Annotation feature leakage')

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
        checks.append(dict(cell_type=c, precision=precision, recall=recall, n_true=int((labels[held] == c).sum()), n_predicted=int((pred == c).sum()), gate_pass=precision >= 0.8 and recall >= 0.6))
    save(out / '272_annotation_holdout_metrics.csv', checks)
    save(out / '272_annotation_genes.csv', pd.DataFrame({'gene_symbol': genes[idx]}))
    writej(out / '272_annotation_scope.json', dict(reference=REFERENCE, held_out_donors=test_donors, method='Multinomial logistic reference label transfer; C=1, balanced classes', probability_threshold=0.7, probability_gap=0.2, excluded_all_fixed_endpoints=True, diagnosis_used_for_annotation=False, target_annotation_author_confirmed=False, target_accuracy_known=False, subtype_assignment='not_attempted', reference_holdout_does_not_prove_target_accuracy=True))
    gate = {r['cell_type']: bool(r['gate_pass']) for r in checks}
    need(gate[astro[0]], 'Astrocyte reference holdout gate failed; no target-label hypothesis testing')
    cpred = centroid_predict(z[held][:, idx], means[:, idx], types)
    centroid_checks = []
    for c in types:
        tp = int(((cpred == c) & (labels[held] == c)).sum())
        precision = tp / max(1, int((cpred == c).sum()))
        recall = tp / max(1, int((labels[held] == c).sum()))
        centroid_checks.append(dict(cell_type=c, precision=precision, recall=recall, gate_pass=precision >= 0.8 and recall >= 0.6))
    save(out / '285_centroid_reference_holdout.csv', centroid_checks)
    model = fit(sample_rows(np.ones(len(obs), dtype=bool)))
    model.audit_centers_ = np.stack([np.asarray(z[labels == c][:, idx].mean(axis=0)).ravel() for c in types])
    model.audit_classes_ = types
    model.audit_gate_ = next((r['gate_pass'] for r in centroid_checks if r['cell_type'] == astro[0]))
    save(out / '285_centroid_reference_centers.csv', pd.DataFrame(model.audit_centers_, index=types, columns=genes[idx]).reset_index(names='reference_class'))
    return (model, genes[idx], astro[0], oligo[0], gate)

def centroid_predict(x, centers, classes):
    norms = np.linalg.norm(centers, axis=1)
    need(np.isfinite(centers).all() and (norms > 0).all(), 'Invalid reference centroids')
    scores = np.asarray(x @ (centers / norms[:, None]).T)
    chosen = np.asarray(classes, dtype=object)[scores.argmax(axis=1)]
    chosen[np.asarray(x.sum(axis=1)).ravel() == 0] = 'Uncertain'
    return chosen

def classify_prob(classes, prob):
    order = np.argsort(prob, axis=1)
    confidence = prob[np.arange(len(prob)), order[:, -1]]
    gap = confidence - prob[np.arange(len(prob)), order[:, -2]]
    pred = np.asarray(classes, dtype=object)[order[:, -1]].copy()
    pred[(confidence < 0.7) | (gap < 0.2)] = 'Uncertain'
    return (pred, confidence, gap)

def number(value):
    value = str(value).strip()
    return float(value) if re.fullmatch('\\d+(?:\\.\\d+)?', value) else np.nan

def normalized_id(value):
    return re.sub('[\\s/_-]', '', str(value)).upper()

def plain_soft(path):
    entities, current = ([], None)
    with Path(path).open('rt', encoding='utf-8-sig') as stream:
        for line in stream:
            line = line.rstrip('\r\n')
            if line.startswith('^') and ' = ' in line:
                kind, value = line[1:].split(' = ', 1)
                current = dict(kind=kind, id=value, fields={})
                entities.append(current)
            elif current is not None and line.startswith('!') and (' = ' in line):
                key, value = line[1:].split(' = ', 1)
                current['fields'].setdefault(key, []).append(value)
    need(any((e['kind'] == 'SAMPLE' for e in entities)), 'Expected GEO sample metadata, not HTML')
    return entities

def target_register(entities):
    rows = []
    for entity in entities:
        if entity['kind'] != 'SAMPLE':
            continue
        c = chars(entity)
        title = first(entity, 'Sample_title')
        fraction = {'parenchyma': 'P', 'vascular': 'V'}.get(c.get('fraction', '').lower())
        need(fraction is not None, f'Unknown fraction {entity['id']}: {c.get('fraction')}')
        need(re.fullmatch(f'{fraction}_\\d{{2}}', title), f'Title/fraction conflict {entity['id']}')
        urls = [u for key, values in entity['fields'].items() if key == 'Sample_supplementary_file' or key.startswith('Sample_supplementary_file_') for u in values]
        urls = [u.replace('ftp://', 'https://', 1) for u in urls if u.endswith('.h5')]
        need(len(urls) == 1, f'Expected one H5 for {entity['id']}')
        need(urllib.parse.urlparse(urls[0]).hostname == 'ftp.ncbi.nlm.nih.gov', 'Unexpected count host')
        name = urllib.parse.unquote(urls[0].rsplit('/', 1)[-1])
        need(name == f'{entity['id']}_{title}_filtered_feature_bc_matrix.h5', 'GSM/title/H5 mismatch')
        need(c.get('disease') in ('AD', 'CT') and c.get('sex') in ('F', 'M'), 'Unexpected diagnosis/sex')
        age = number(c.get('age', ''))
        need(np.isfinite(age) and 18 <= age <= 110, 'Missing or unparseable age')
        donor = c.get('donor id', '').strip()
        need(bool(donor), 'Explicit donor ID required')
        need('frontal cortex' in c.get('tissue', '').lower(), 'Unexpected tissue')
        processing = ' '.join(entity['fields'].get('Sample_data_processing', []))
        need('CellRanger' in processing and 'filtered' in processing and ('HDF5' in processing), 'Source count semantics changed')
        rows.append(dict(gsm=entity['id'], title=title, fraction=fraction, donor=donor, case=int(c['disease'] == 'AD'), age=age, sex=c['sex'], pmi=number(c.get('pmi', '')), ph_cerebellum=number(c.get('ph cerebellum', '')), rin=number(c.get('rin', '')), braak=c.get('braak tangle', ''), apoe=c.get('apoe', ''), tissue=c['tissue'], count_url=urls[0], source_name=first(entity, 'Sample_source_name_ch1'), processing=processing, original_characteristics=json.dumps(c, sort_keys=True)))
    m = pd.DataFrame(rows)
    need(len(m) == 40 and m.gsm.is_unique and m.title.is_unique, 'Expected exactly 40 unique GEO libraries')
    need(m.donor.nunique() == 20 and (not m.duplicated(['donor', 'fraction']).any()), 'Expected 20 donor pairs')
    for donor, sub in m.groupby('donor'):
        need(set(sub.fraction) == {'P', 'V'}, f'Incomplete fractions: {donor}')
        for col in ['case', 'age', 'sex', 'pmi', 'ph_cerebellum', 'rin', 'braak', 'apoe']:
            need(sub[col].fillna('__missing__').nunique() == 1, f'Conflicting donor field: {donor}/{col}')
    need(m.drop_duplicates('donor').case.value_counts().to_dict() == {0: 10, 1: 10}, 'Expected 10 AD / 10 CT')
    return m.sort_values(['fraction', 'title']).reset_index(drop=True)

def read_counts(path):
    with h5py.File(path, 'r') as f:
        need('matrix' in f, 'Not a Cell Ranger v3 H5')
        n = f['matrix']
        features = n['features']
        symbols, ids = (decoded(features['name'][:]), decoded(features['id'][:]))
        types, barcodes = (decoded(features['feature_type'][:]), decoded(n['barcodes'][:]))
        shape = tuple(n['shape'][:])
        data = n['data'][:]
        need(shape == (len(symbols), len(barcodes)), 'H5 dimensions mismatch')
        need(len(set(ids)) == len(ids) and len(set(barcodes)) == len(barcodes), 'Duplicate feature IDs/barcodes')
        need(np.isfinite(data).all() and (data >= 0).all() and (data == np.floor(data)).all(), 'Not integer UMI counts')
        need((data < 2 ** 53).all(), 'Counts too large to convert exactly')
        x = sp.csc_matrix((data.astype(np.int64), n['indices'][:], n['indptr'][:]), shape=shape)
        x.check_format(full_check=True)
    use = types == 'Gene Expression'
    need(use.sum() >= 1000, 'Too few gene features')
    symbols = symbols[use]
    ids = ids[use]
    need(all((str(s).strip() for s in symbols)), 'Empty gene symbol')
    unique, inv = np.unique(symbols, return_inverse=True)
    mapper = sp.csr_matrix((np.ones(len(inv), dtype=np.int64), (np.arange(len(inv)), inv)), shape=(len(inv), len(unique)))
    y = (x[use].T @ mapper).tocsr()
    y.eliminate_zeros()
    return (y, unique, barcodes, pd.DataFrame(dict(feature_id=ids, gene_symbol=symbols)))

class GeneIdentity:
    """Approved IDs resolve identities; ambiguous historical tokens are never guessed."""

    def __init__(self, frame):
        required = {'hgnc_id', 'symbol', 'status', 'prev_symbol', 'alias_symbol', 'ensembl_gene_id'}
        need(required <= set(frame), 'HGNC schema incomplete')
        self.frame = frame.fillna('').astype(str)
        self.frame = self.frame[self.frame.status.eq('Approved')].copy()
        need(self.frame.hgnc_id.is_unique and self.frame.symbol.is_unique, 'Duplicate approved HGNC IDs/symbols')
        self.rows = self.frame.set_index('hgnc_id').to_dict('index')
        self.approved, self.historical, self.ensembl, self.all_tokens = ({}, {}, {}, {})
        for row in self.frame.itertuples(index=False):
            self.approved.setdefault(row.symbol.upper(), set()).add(row.hgnc_id)
            for field in ('prev_symbol', 'alias_symbol'):
                for token in getattr(row, field).split('|'):
                    if token:
                        self.historical.setdefault(token.upper(), set()).add(row.hgnc_id)
            for token in row.ensembl_gene_id.split('|'):
                if token:
                    self.ensembl.setdefault(token.split('.')[0], set()).add(row.hgnc_id)
            tokens = {row.symbol.upper()}
            for field in ('prev_symbol', 'alias_symbol'):
                tokens.update((t.upper() for t in getattr(row, field).split('|') if t))
            self.all_tokens[row.hgnc_id] = tokens

    def candidates(self, symbol):
        key = str(symbol).strip().upper()
        if key in self.approved:
            return (set(self.approved[key]), 'approved_casefold')
        return (set(self.historical.get(key, set())), 'historical_symbol')

    def exclusion_tokens(self, symbols):
        tokens = {str(s).strip().upper() for s in symbols}
        for symbol in list(tokens):
            ids = self.approved.get(symbol, set()) | self.historical.get(symbol, set())
            for identity in ids:
                tokens.update(self.all_tokens[identity])
        return tokens

    def target_register(self, features):
        rows = []
        for symbol, group in features.groupby('gene_symbol', sort=True):
            symids, _ = self.candidates(symbol)
            ensgs = sorted(set((str(x).split('.')[0] for x in group.feature_id)))
            known = [self.ensembl.get(x, set()) for x in ensgs]
            stable = set().union(*known)
            identity = ''
            reason = 'unresolved_target_identity'
            if len(stable) == 1 and all((len(x) == 1 for x in known)):
                identity = next(iter(stable))
                reason = 'ensembl_id'
                if symids and identity not in symids:
                    identity = ''
                    reason = 'ensembl_symbol_conflict'
            elif stable:
                reason = 'ambiguous_or_partially_resolved_feature_ids'
            elif len(symids) == 1:
                identity = next(iter(symids))
                reason = 'unique_hgnc_symbol_without_current_ensembl_link'
            elif len(symids) > 1:
                reason = 'ambiguous_target_symbol'
            rows.append(dict(target_symbol=symbol, feature_ids='|'.join(ensgs), hgnc_id=identity, approved_symbol=self.rows[identity]['symbol'] if identity else '', target_identity_status=reason))
        return pd.DataFrame(rows)

    def map_endpoints(self, fixed, targets):
        rows = []
        for original in fixed.gene_symbol:
            ids, method = self.candidates(original)
            identity = next(iter(ids)) if len(ids) == 1 else ''
            found = targets[targets.hgnc_id.eq(identity)] if identity else targets.iloc[:0]
            status = 'mapped' if len(found) == 1 else 'target_identity_absent' if identity and found.empty else 'multiple_target_rows' if identity else 'ambiguous_endpoint' if ids else 'unresolved_endpoint'
            rows.append(dict(gene_symbol=original, hgnc_id=identity, approved_symbol=self.rows[identity]['symbol'] if identity else '', analysis_symbol=str(found.iloc[0].target_symbol) if status == 'mapped' else '', mapping_status=status, mapping_method=method, candidate_hgnc_ids='|'.join(sorted(ids))))
        result = pd.DataFrame(rows)
        collision = result.hgnc_id.ne('') & result.hgnc_id.duplicated(keep=False)
        result.loc[collision, 'mapping_status'] = 'duplicate_fixed_identity'
        result.loc[collision, 'analysis_symbol'] = ''
        need(result.loc[result.mapping_status.eq('mapped'), 'analysis_symbol'].is_unique, 'Target mapping collision')
        return result

def obtain_identity(dl, out, root, explicit=None):
    if explicit:
        path = Path(explicit).resolve()
        receipt = path.with_name(path.name + '.receipt.json')
        need(path.is_file() and receipt.is_file(), 'Local HGNC requires its official URL/hash .receipt.json')
        r = readj(receipt)
        need(r.get('url') in HGNC_URLS and r.get('sha256') == digest(path), 'HGNC receipt mismatch')
        dl.receipts.append(dict(r, local_path=str(path)))
    else:
        path = None
        errors = []
        for i, url in enumerate(HGNC_URLS):
            try:
                path = dl.get(url, 'HGNC_complete_' + str(i) + '.tsv')
                break
            except (Blocked, OSError) as e:
                errors.append(str(e))
        need(path is not None, 'Official HGNC unavailable; ' + str(errors))
    table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    need(len(table) >= 30000, 'HGNC snapshot is incomplete/unexpected')
    identity = GeneIdentity(table)
    sha = digest(path)
    lock = dict(sha256=sha, rows=len(table), source_url=next((r['url'] for r in reversed(dl.receipts) if r.get('sha256') == sha)))
    lockpath = root / '284_gene_identity_lock.json'
    if lockpath.exists():
        need(readj(lockpath) == lock, 'HGNC snapshot changed; use a new output name')
    else:
        writej(lockpath, lock)
    shutil.copyfile(path, out / '284_HGNC_complete.tsv')
    writej(out / '284_gene_identity_source.json', lock)
    return identity

def qc_counts(x, genes):
    umi = np.asarray(x.sum(axis=1)).ravel()
    detected = np.diff(x.indptr)
    mito = np.asarray(x[:, np.char.startswith(genes.astype(str), 'MT-')].sum(axis=1)).ravel()
    frac = np.divide(mito, umi, out=np.zeros(len(umi), float), where=umi > 0)
    keep = (umi >= POLICY['min_umi']) & (umi <= POLICY['max_umi']) & (detected >= POLICY['min_features']) & (detected <= POLICY['max_features']) & (frac <= POLICY['max_mito'])
    return (keep, umi, detected, frac)

def annotate_and_aggregate(samples, paths, dl, fixed, out, identity):
    print('285: read author-labeled reference and train endpoint-excluding classifier', flush=True)
    readme = dl.get(geo(REFERENCE) + 'suppl/GSE243292_processed_data_readme.txt')
    shutil.copyfile(readme, out / '285_reference_processed_data_readme.txt')
    need('Count matrix' in readme.read_text(encoding='utf-8-sig'), 'Reference source count declaration changed')
    first_x, first_genes, _, first_mapping = read_counts(paths[samples.iloc[0].gsm])
    del first_x
    targets = identity.target_register(first_mapping)
    save(out / '285_target_gene_identity.csv', targets)
    mapping = identity.map_endpoints(fixed, targets)
    mapped_fixed = fixed.merge(mapping, on='gene_symbol', how='left', validate='one_to_one')
    save(out / '285_fixed_endpoint_identity.csv', mapped_fixed)
    save(out / '285_identity_coverage_summary.csv', mapped_fixed.groupby('mapping_status', sort=True).size().reset_index(name='endpoint_slots'))
    rx, rg, obs = get_reference(dl)
    tokens = identity.exclusion_tokens(list(fixed.gene_symbol) + sum(MARKERS.values(), []))
    forbidden = {str(g) for g in set(rg) | set(first_genes) if str(g).upper() in tokens}
    forbidden.update(mapped_fixed.analysis_symbol[mapped_fixed.analysis_symbol.ne('')])
    model, features, astro, oligo, gates = classifier(rx, rg, obs, forbidden, out, seed=284, allowed=set(first_genes))
    need(not any((str(g).upper() in tokens or g in forbidden for g in features)), 'Fixed endpoint/heldout marker identity leaked into annotation features')
    save(out / '285_annotation_identity_exclusions.csv', [dict(gene_symbol=g, reason='fixed_endpoint_or_heldout_marker_identity') for g in sorted(forbidden)])
    weights = pd.DataFrame(model.coef_, index=model.classes_, columns=features)
    weights.index.name = 'reference_class'
    weights.to_csv(out / '285_classifier_coefficients.csv')
    save(out / '285_classifier_intercepts.csv', pd.DataFrame(dict(reference_class=model.classes_, intercept=model.intercept_)))
    del rx
    writej(out / '285_annotation_provenance.json', dict(reference=REFERENCE, reference_original_label_field='final_celltype', label_type='new_reference_transferred_inferred_labels_not_target_author_labels', target_author_barcode_annotation_obtained=False, target_accuracy_known=False, labels_linked_by=['GSM', 'original_barcode'], forbidden=sorted(forbidden), primary_fraction='P', secondary_fraction='V', never_pool_fractions=True, primary_label_mode='fresh_fit', ood_limit='Reference broad classes do not cover every vascular subtype; confidence is not target validation', classifier_classes=list(map(str, model.classes_)), reference_holdout_gates=gates))
    cols, units, profiles, gene_master = ([], [], [], None)
    consensus_cols, consensus_units, agreement = ([], [], [])
    for row in samples.itertuples(index=False):
        print(f'285: annotate {row.gsm} / {row.title}', flush=True)
        x, genes, barcodes, mapping = read_counts(paths[row.gsm])
        a = mapping.sort_values(['feature_id', 'gene_symbol']).reset_index(drop=True)
        b = first_mapping.sort_values(['feature_id', 'gene_symbol']).reset_index(drop=True)
        need(a.equals(b), 'Feature-ID/symbol correspondence differs between target libraries')
        if gene_master is None:
            gene_master = genes
            save(out / '285_feature_symbol_mapping.csv', mapping)
        need(np.array_equal(gene_master, genes), 'Gene universe differs between target libraries; no silent zero-fill')
        idx = pd.Index(genes).get_indexer(features)
        need((idx >= 0).all(), 'Missing annotation features')
        keep, umi, detected, mito = qc_counts(x, genes)
        labels = np.full(len(barcodes), 'QC_excluded', object)
        confidence, gap = (np.full(len(barcodes), np.nan), np.full(len(barcodes), np.nan))
        ids = np.flatnonzero(keep)
        centroid_labels = np.full(len(barcodes), 'QC_excluded', object)
        if len(ids):
            y = lognorm(x[ids])
            probabilities = model.predict_proba(y[:, idx])
            pred, conf, margin = classify_prob(model.classes_, probabilities)
            labels[ids], confidence[ids], gap[ids] = (pred, conf, margin)
            centroid_labels[ids] = centroid_predict(y[:, idx], model.audit_centers_, model.audit_classes_)
            del y
        marker_idx = {k: np.flatnonzero(np.isin(genes, values)) for k, values in MARKERS.items()}
        marker_n = {k: np.asarray((x[:, ix] > 0).sum(axis=1)).ravel() for k, ix in marker_idx.items()}
        conflict = (labels == astro) & (marker_n['astro'] == 0) & ((marker_n['vascular'] >= 2) | (marker_n['oligo'] >= 2))
        final = labels.copy()
        final[conflict] = 'Uncertain_marker_conflict'
        table = pd.DataFrame(dict(gsm=row.gsm, barcode=barcodes, donor=row.donor, fraction=row.fraction, qc_pass=keep, total_umi=umi, detected_genes=detected, mito_fraction=mito, inferred_label=labels, label=final, confidence=confidence, probability_gap=gap, heldout_marker_conflict=conflict))
        refit_changes = 0
        table['centroid_label'] = centroid_labels
        table['consensus_astro'] = (final == astro) & (centroid_labels == astro)
        save(out / f'285_{row.gsm}_cell_annotations.csv.gz', table)
        consensus = table.consensus_astro.to_numpy()
        summed = np.asarray(x[consensus].sum(axis=0)).ravel().astype(np.int64)
        consensus_cols.append(summed)
        consensus_units.append(dict(unit_id=f'{row.gsm}__Ast_consensus', gsm=row.gsm, donor=row.donor, fraction=row.fraction, cell_type='Ast_consensus', n_cells=int(consensus.sum()), total_umi=int(summed.sum()), **{c: getattr(row, c) for c in ['case', 'age', 'sex', 'pmi', 'ph_cerebellum', 'rin']}))
        agreement.append(dict(gsm=row.gsm, donor=row.donor, fraction=row.fraction, case=row.case, primary_astro=int((final == astro).sum()), centroid_astro=int((centroid_labels == astro).sum()), consensus_astro=int(consensus.sum()), reference_centroid_gate=bool(model.audit_gate_), diagnostic_refit_changed_labels=refit_changes, primary_retention=float(consensus.sum() / max(1, (final == astro).sum())), boundary='same_reference_algorithm_agreement_not_target_accuracy'))
        for kind in sorted(set(final)):
            if kind == 'QC_excluded':
                continue
            mask = final == kind
            sums = np.asarray(x[mask].sum(axis=0)).ravel().astype(np.int64)
            uid = f'{row.gsm}__{kind}'
            cols.append(sums)
            units.append(dict(unit_id=uid, gsm=row.gsm, donor=row.donor, fraction=row.fraction, cell_type=kind, n_cells=int(mask.sum()), total_umi=int(sums.sum()), **{c: getattr(row, c) for c in ['case', 'age', 'sex', 'pmi', 'ph_cerebellum', 'rin']}))
        ast = final == astro
        other = keep & ~np.isin(final, [astro, 'Uncertain', 'Uncertain_marker_conflict'])
        for kind, ix in marker_idx.items():
            for j in ix:

                def mean_normalized(mask):
                    if not mask.any():
                        return np.nan
                    return float(np.log1p(10000 * x[mask, j].toarray().ravel() / umi[mask]).mean())
                av, ov = (mean_normalized(ast), mean_normalized(other))
                profiles.append(dict(gsm=row.gsm, fraction=row.fraction, marker=genes[j], marker_set=kind, ast_mean_log1p=av, other_mean_log1p=ov, difference=av - ov, ast_cells=int(ast.sum()), other_cells=int(other.sum()), conflict_fraction=float(conflict.sum() / max(1, (labels == astro).sum()))))
    metadata = pd.DataFrame(units)
    need(len(cols) > 0, 'No target cells retained by the locked QC; no expression models')
    counts = pd.DataFrame(np.stack(cols, axis=1), index=gene_master, columns=metadata.unit_id)
    counts.index.name = 'gene_symbol'
    counts.to_csv(out / '286_all_class_pseudobulk_counts.csv.gz')
    cc = pd.DataFrame(np.stack(consensus_cols, axis=1), index=gene_master, columns=[u['unit_id'] for u in consensus_units])
    cc.index.name = 'gene_symbol'
    cc.to_csv(out / '286_consensus_pseudobulk_counts.csv.gz')
    save(out / '286_consensus_units.csv', consensus_units)
    save(out / '285_target_annotation_agreement.csv', agreement)
    save(out / '286_analysis_units.csv', metadata)
    p = save(out / '285_heldout_marker_diagnostics.csv', profiles)
    checks = []
    for fraction in ['P', 'V']:
        sub = p[(p.fraction == fraction) & (p.ast_cells >= POLICY['min_cells']) & (p.other_cells >= 50)]
        diffs = sub[sub.marker_set == 'astro'].groupby('marker').difference.median()
        positives = int((diffs >= POLICY['marker_log_difference']).sum())
        contamination = sub.groupby('gsm').conflict_fraction.first()
        passed = len(contamination) >= 2 * POLICY['min_group'] and positives >= POLICY['min_marker_genes'] and (not contamination.empty) and (contamination.max() <= POLICY['max_contaminant_fraction'])
        checks.append(dict(fraction=fraction, positive_heldout_markers=positives, assessed_libraries=len(contamination), passed=bool(passed), interpretation='Operational target plausibility check, not gold-standard accuracy'))
    save(out / '285_annotation_gate.csv', checks)
    return (counts, metadata, astro, {r['fraction']: r['passed'] for r in checks}, mapped_fixed)

def age_review(m):
    a = m[m.case == 1].age.to_numpy(float)
    b = m[m.case == 0].age.to_numpy(float)
    if min(len(a), len(b)) < 2:
        return dict(pass_gate=False, reason='insufficient_groups')
    lo, hi = (max(a.min(), b.min()), min(a.max(), b.max()))
    retention = [float(((x >= lo) & (x <= hi)).mean()) if hi >= lo else 0.0 for x in [a, b]]
    sd = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)
    smd = float((a.mean() - b.mean()) / sd) if sd > 0 else 0.0 if a.mean() == b.mean() else None
    bins = np.arange(0, 126, 5)
    overlap = float(np.minimum(np.histogram(a, bins)[0] / len(a), np.histogram(b, bins)[0] / len(b)).sum())
    ok = min(retention) >= POLICY['min_age_range_retention'] and smd is not None and (abs(smd) <= POLICY['max_abs_age_smd']) and (overlap >= POLICY['min_age_bin_overlap'])
    return dict(pass_gate=bool(ok), n_AD=len(a), n_CT=len(b), overlap_low=float(lo), overlap_high=float(hi), AD_retention=retention[0], CT_retention=retention[1], SMD=smd, age_bin_overlap=overlap, reason='operational_rule_not_proof_of_conditional_overlap')
