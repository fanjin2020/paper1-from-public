"""Public-input supplementary scientific kernels."""
from pathlib import Path
import csv, gzip, io, json, os, re, shutil, subprocess
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
from pipeline.io import require as need, sha, dump as writej
from .common import csvout, readj, reserved_bh as bh

js = writej


SOURCES = {
    'GSE1297_series_matrix.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1297/matrix/GSE1297_series_matrix.txt.gz',
    'GPL96.annot.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz',
}


def table(text, start, end):
    need(text.count(start) == 1 and text.count(end) == 1, 'Unexpected GEO table markers')
    body = text.split(start)[1].split(end)[0].strip()
    return pd.read_csv(io.StringIO(body), sep='\t', dtype=str, keep_default_na=False)


def prepare(matrix_file, annot_file, out):
    text = gzip.decompress(matrix_file.read_bytes()).decode('utf-8-sig')
    rows = list(csv.reader(io.StringIO(text.split('!series_matrix_table_begin')[0]), delimiter='\t'))
    fields = {}
    for row in rows:
        if row:
            fields.setdefault(row[0], []).append(row[1:])
    need(fields.get('!Series_geo_accession') == [['GSE1297']], 'Wrong GEO series')
    ids = fields['!Sample_geo_accession'][0]
    need(len(ids) == 31 and len(set(ids)) == 31, 'Expected 31 unique arrays')
    need(set(fields['!Sample_source_name_ch1'][0]) == {'hippocampal CA1 tissue'}, 'Unexpected tissue')
    need(set(fields['!Sample_platform_id'][0]) == {'GPL96'}, 'Unexpected platform')
    records = [dict(sample_id=g) for g in ids]
    for values in fields['!Sample_characteristics_ch1']:
        need(len(values) == len(ids), 'Characteristic length mismatch')
        for rec, value in zip(records, values):
            key, sep, val = value.partition(':')
            need(bool(sep), f'Unrecognized characteristic: {value}')
            key = key.strip().lower()
            need(key not in rec, f'Duplicate field: {key}')
            rec[key] = val.strip()
    m = pd.DataFrame(records).rename(columns={'ref #': 'source_subject_label'})
    need({'mmse', 'braak', 'nft', 'age', 'sex', 'pmi', 'source_subject_label'} <= set(m), 'Missing clinical fields')
    need(not m.source_subject_label.duplicated().any(), 'Repeated source subject label; do not treat repeats as donors')
    for col in ['mmse', 'braak', 'nft', 'age', 'pmi']:
        m[col] = pd.to_numeric(m[col], errors='raise')
        need(np.isfinite(m[col]).all(), f'Incomplete {col}; no imputation permitted')
    need(m.mmse.between(0, 30).all() and m.braak.between(0, 6).all(), 'Score outside defined scale')
    need(m.age.between(18, 115).all() and m.nft.ge(0).all() and m.pmi.ge(0).all(), 'Invalid covariate')
    need(m.sex.isin(['M', 'F']).all(), 'Unknown sex encoding')
    m['score_semantics'] = 'published_MMSE_may_include_author_time_adjustment'
    m.to_csv(out / '230_all_samples.csv', index=False)
    pd.DataFrame([{'field': k, 'value': '\t'.join(v)} for k, values in fields.items() for v in values]).to_csv(out / '230_source_declarations.csv', index=False)
    t = table(text, '!series_matrix_table_begin', '!series_matrix_table_end').set_index('ID_REF')
    need(t.columns.tolist() == ids and t.index.is_unique, 'Matrix identifiers/order mismatch')
    y = t.apply(pd.to_numeric, errors='raise')
    atext = gzip.decompress(annot_file.read_bytes()).decode('utf-8-sig')
    annot = table(atext, '!platform_table_begin', '!platform_table_end')
    need(annot.ID.is_unique and {'ID', 'Gene symbol'} <= set(annot), 'Invalid platform annotation')
    mapping = annot.set_index('ID')['Gene symbol'].reindex(y.index).fillna('')
    valid_symbol = mapping.str.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*')
    valid_values = pd.Series(np.isfinite(y.to_numpy()).all(axis=1), index=y.index) & y.gt(0).all(axis=1)
    decision = pd.DataFrame({'probe_id': y.index, 'raw_gene_symbol': mapping.values,
                             'single_symbol': valid_symbol.values, 'positive_complete': valid_values.values})
    decision['included'] = decision.single_symbol & decision.positive_complete
    decision.to_csv(out / '230_probe_background_decisions.csv', index=False)
    keep = valid_symbol & valid_values
    need(int(keep.sum()) >= 1000, 'Insufficient eligible probes')
    # Published MAS5 signal intensities: log2, then unweighted mean across unique-symbol probes.
    # No clinical-outcome-dependent probe selection or additional quantile normalization.
    logy = np.log2(y.loc[keep])
    logy.index = mapping.loc[keep].values
    gene = logy.groupby(level=0, sort=True).mean()
    gene.index.name = 'gene_symbol'
    variable = gene.var(axis=1).gt(1e-12)
    pd.DataFrame({'gene_symbol': gene.index, 'included_nonconstant': variable.values}).to_csv(out / '230_gene_background_decisions.csv', index=False)
    gene.loc[variable].to_csv(out / '230_expression_log2.csv.gz', compression={'method': 'gzip', 'compresslevel': 1})
    pd.Series(mapping.loc[keep].value_counts(), name='n_probes').rename_axis('gene_symbol').to_csv(out / '230_gene_probe_counts.csv')
    js(out / '230_preparation.json', {'n_samples': len(m), 'input_probes': len(y), 'retained_probes': int(keep.sum()),
       'eligible_genes': int(variable.sum()), 'scale': 'log2 deposited positive MAS5 intensities; mean log2 per gene',
       'preprocessing_limit': 'Published processing and cohort selection cannot be undone; not a raw-CEL RMA reanalysis'})
