"""Public-input supplementary scientific kernels."""
from pathlib import Path
import csv, gzip, io, json, os, re, shutil, subprocess
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
from pipeline.io import require as need, sha, dump as writej
from .common import csvout, readj, reserved_bh as bh
import contextlib, importlib.metadata, traceback
from . import interaction_math as U
from .common import save, Blocked
require = need
csv = pd.read_csv


SEEDS = (291, 292)


COMPARATORS = ('End', 'Ex', 'In', 'Mic', 'Oli', 'Opc')


FOCAL = ('GFAP', 'GSDMD', 'IGFBP5')


class NotTestable(RuntimeError):
    pass


def boolean(values, label):
    require(values.notna().all() and values.isin([True, False]).all(), 'Invalid Boolean column: ' + label)
    return values.to_numpy(dtype=bool)


def score_alignment(observed, ann, keep, threshold):
    require(observed.barcode.is_unique and observed.barcode.notna().all(), 'Duplicate/missing score barcode')
    require(ann.barcode.is_unique and ann.barcode.notna().all(), 'Duplicate/missing annotation barcode')
    wanted = ann.loc[keep, 'barcode'].tolist()
    require(set(observed.barcode) == set(wanted), 'Score barcodes do not equal all QC-passing barcodes')
    obs = observed.set_index('barcode').loc[wanted]
    require(np.isfinite(obs.score).all() and 0 < threshold < 1, 'Invalid score or threshold')
    flags = boolean(obs.predicted_doublet, 'predicted_doublet')
    require(np.array_equal(flags, obs.score.to_numpy() > threshold), 'Recorded calls differ from original threshold')
    retained = np.zeros(len(ann), bool)
    retained[np.flatnonzero(keep)] = ~flags
    return retained, obs.reset_index()


def fixed_grid(results):
    variants = ('unfiltered_matched', 'filtered_matched')
    grid = pd.MultiIndex.from_product([variants, COMPARATORS, results['endpoints'].gene_symbol],
        names=['variant', 'comparator', 'gene_symbol']).to_frame(index=False)
    frames = results['frames']
    if frames:
        z = pd.concat(frames, ignore_index=True)
        require(not z.duplicated(['variant', 'comparator', 'gene_symbol']).any(), 'Duplicate fixed result')
        grid = grid.merge(z, on=['variant', 'comparator', 'gene_symbol'], how='left', validate='one_to_one')
    for field in ['P_Value', 'logFC', 'CI_low', 'CI_high', 'SE_moderated', 'df_total', 'n_donors']:
        if field not in grid:
            grid[field] = np.nan
    grid['testable'] = grid.P_Value.notna()
    grid['q_1422_within_variant'] = grid.groupby('variant', sort=False).P_Value.transform(U.reserved_bh)
    grid['q_2844_all_variants_supplement'] = U.reserved_bh(grid.P_Value)
    n = len(results['endpoints']) * len(COMPARATORS)
    grid['primary_reporting_family'] = f'filtered_matched:{n};unfiltered_matched:{n}_sensitivity;{2*n}_supplement'
    grid['within_variant_family_size'] = n
    grid['joint_family_size'] = 2 * n
    grid['scope'] = 'post_selection_sensitivity_not_independent_validation'
    require(len(grid) == 2 * n, 'Incomplete planned family')
    return grid


def common_background_engine(comparator):
    script, arms = U.boundary_engine('interaction_second_rho', comparator)
    anchor = "ko<-filterByExpr(DGEList(Y[,io,drop=FALSE]),group=m$case[io]);keep<-ka&ko"
    require(script.count(anchor) == 1, 'Reviewed R template anchor changed')
    replacement = """ko<-filterByExpr(DGEList(Y[,io,drop=FALSE]),group=m$case[io])
 Other<-readcount('other_counts.csv.gz')
 need(identical(dimnames(Y),dimnames(Other)), 'Paired filter matrices differ in order')
 koa<-filterByExpr(DGEList(Other[,ia,drop=FALSE]),group=m$case[ia])
 koo<-filterByExpr(DGEList(Other[,io,drop=FALSE]),group=m$case[io])
 keep<-ka&ko&koa&koo
 write.csv(data.frame(gene_symbol=rownames(Y),current_ast=ka,current_other=ko,
 paired_ast=koa,paired_other=koo,common=keep),file.path(d,'filter_comparison.csv'),row.names=FALSE)"""
    script = script.replace(anchor, replacement)
    gate = "need(sum(keep)>=500,'Fewer than 500 interaction genes')"
    require(script.count(gate) == 1, 'Background gate anchor changed')
    script = script.replace(gate, """write.csv(data.frame(common_genes=sum(keep),minimum=500,passed=sum(keep)>=500),
 file.path(d,'common_background_gate.csv'),row.names=FALSE)
 if(sum(keep)<500)quit(save='no',status=42)""")
    return script, arms


def execute_pair_model(d, comparator, fixed, rscript):
    script, arms = common_background_engine(comparator)
    (d / 'engine.R').write_text(script, encoding='utf-8')
    with (d / 'engine.log').open('w', encoding='utf-8') as log:
        rc = subprocess.run([rscript, '--vanilla', str(d / 'engine.R'), str(d), 'interaction'],
                            stdout=log, stderr=subprocess.STDOUT).returncode
    if rc == 42:
        gate = csv(d / 'common_background_gate.csv')
        require(len(gate) == 1 and int(gate.minimum.iloc[0]) == 500 and
                int(gate.common_genes.iloc[0]) < 500 and not bool(gate.passed.iloc[0]), 'Invalid R not-testable gate')
        raise NotTestable('Fewer than 500 genes in both classes and count variants; no background relaxation')
    require(rc == 0, 'R model failed: ' + str(d / 'engine.log'))
    return U.verify_sensitivity_outputs(d, 'interaction_second_rho', fixed, expected_arms=arms)


def paired_setup(primary, units, comparator):
    p = units[(units.fraction == 'P') & units.cell_type.isin(['Ast', comparator])].copy()
    require(not p.duplicated(['donor', 'cell_type']).any(), 'Duplicate donor/class units')
    eligible = p[p.n_cells >= U.POLICY['min_cells']]
    ids = set(eligible.loc[eligible.cell_type == 'Ast', 'donor']) & set(eligible.loc[eligible.cell_type == comparator, 'donor'])
    m = primary[primary.donor.isin(ids)].copy().reset_index(drop=True)
    B = U.matching_gate(m)
    indexed = p.set_index(['donor', 'cell_type'])
    a = indexed.loc[[(donor, 'Ast') for donor in m.donor]].reset_index()
    o = indexed.loc[[(donor, comparator) for donor in m.donor]].reset_index()
    for table in [a, o]:
        for field in ['case', 'age', 'sex', 'pmi', 'ph_cerebellum']:
            require(np.array_equal(table[field].to_numpy(), m[field].to_numpy()), 'Donor covariate mismatch: ' + field)
    pairs = pd.concat([o, a], ignore_index=True)
    bx = np.tile(B.to_numpy(), (2, 1))
    ast = np.repeat([0., 1.], len(m))
    X = pd.DataFrame(np.column_stack([bx, bx * ast[:, None]]), index=pairs.unit_id,
                     columns=list(B.columns) + ['ast_' + c for c in B.columns])
    U.need(np.linalg.matrix_rank(X) == X.shape[1] and len(X) - X.shape[1] >= U.POLICY['min_df'], 'Paired design not estimable')
    X.index.name = 'unit_id'
    return pairs, X


def seed_diagnostic(x, keep, original, sample, d):
    import scrublet
    from scipy.stats import spearmanr
    from threadpoolctl import threadpool_limits
    d.mkdir(parents=True)
    y = x[keep].tocsc()
    y = y[:, np.asarray(y.sum(axis=0)).ravel() > 0].tocsr()
    rows = []
    for seed in SEEDS:
        prefix = d / str(seed)
        try:
            with prefix.with_suffix('.log').open('w', encoding='utf-8') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log), threadpool_limits(limits=1):
                sc = scrublet.Scrublet(y, expected_doublet_rate=.06, random_state=seed)
                scores, calls = sc.scrub_doublets(use_approx_neighbors=False, n_prin_comps=30, verbose=True)
            scores = np.asarray(scores, float)
            sim = np.asarray(sc.doublet_scores_sim_, float)
            require(scores.shape == (len(original),) and np.isfinite(scores).all() and np.isfinite(sim).all(), 'Invalid repeat score output')
            threshold = getattr(sc, 'threshold_', None)
            usable = calls is not None and threshold is not None and np.isfinite(threshold)
            row = dict(gsm=sample.gsm, fraction=sample.fraction, donor=sample.donor, seed=seed,
                       status='completed_diagnostic' if usable else 'automatic_threshold_unavailable',
                       threshold=float(threshold) if usable else None,
                       expected_rate_assumption=.06, package_version=importlib.metadata.version('scrublet'))
            rho = spearmanr(scores, original.score.to_numpy()).statistic
            row['score_spearman_vs_original'] = float(rho) if np.isfinite(rho) else None
            table = pd.DataFrame({'barcode': original.barcode, 'score': scores})
            if usable:
                calls = np.asarray(calls)
                require(calls.dtype.kind == 'b' and np.array_equal(calls, scores > threshold), 'Repeat threshold/calls mismatch')
                prior = boolean(original.predicted_doublet, 'original predicted_doublet')
                union = int((calls | prior).sum())
                row.update(predicted_fraction=float(calls.mean()), simulated_detectable_fraction=float((sim > threshold).mean()),
                           call_agreement=float((calls == prior).mean()),
                           doublet_jaccard=float((calls & prior).sum() / union) if union else None)
                table['predicted_doublet'] = calls
            save(prefix.with_name(prefix.name + '_observed.csv.gz'), table)
            save(prefix.with_name(prefix.name + '_simulated.csv.gz'), pd.DataFrame({'score': sim}))
            rows.append(row)
        except Exception as e:
            rows.append(dict(gsm=sample.gsm, seed=seed, status='technical_error', reason=str(e)))
            prefix.with_suffix('.error.txt').write_text(traceback.format_exc(), encoding='utf-8')
        save(d / 'summary.csv', rows)
    return rows


def models(source, out, original, filtered, units, rscript):
    ep = csv(source.get('285_fixed_endpoint_identity.csv'), dtype=str, keep_default_na=False)
    require(len(ep) > 0 and ep.gene_symbol.is_unique, 'Invalid current discovery endpoint family')
    primary = csv(source.get('288_P_models/metadata.csv'), dtype={'donor': str, 'unit_id': str})
    hold = csv(source.get('272_annotation_holdout_metrics.csv')).set_index('cell_type')
    require(hold.index.is_unique, 'Duplicate reference class gate')
    frames, statuses = [], []
    for comparator in COMPARATORS:
        print('298 paired filtered interaction:', comparator, flush=True)
        parent = out / ('298_Ast_vs_' + comparator); parent.mkdir()
        try:
            for label in ['Ast', comparator]:
                U.need(label in hold.index and bool(hold.loc[label, 'gate_pass'] == True), 'Reference holdout gate not passed: ' + label)
            a = units[(units.fraction == 'P') & (units.cell_type == 'Ast') & (units.n_cells >= U.POLICY['min_cells'])]
            b = units[(units.fraction == 'P') & (units.cell_type == comparator) & (units.n_cells >= U.POLICY['min_cells'])]
            ids = set(a.donor) & set(b.donor)
            save(parent / 'donor_flow.csv', primary[['donor', 'case', 'unit_id']].assign(included=primary.donor.isin(ids)))
            pairs, X = paired_setup(primary, units, comparator)
            require(set(pairs.unit_id).issubset(original.columns) and set(pairs.unit_id).issubset(filtered.columns), 'Missing paired matrix columns')
        except U.Blocked as e:
            statuses.append(dict(comparator=comparator, status='not_testable', reason=str(e)))
            save(out / '298_model_status.csv', statuses)
            continue
        pair_frames = []
        try:
            for variant, now, other in [('unfiltered_matched', original, filtered), ('filtered_matched', filtered, original)]:
                d = parent / variant; d.mkdir()
                save(d / 'metadata.csv', pairs); X.to_csv(d / 'design.csv'); save(d / 'endpoints.csv', ep)
                now.loc[:, pairs.unit_id].to_csv(d / 'counts.csv.gz')
                other.loc[:, pairs.unit_id].to_csv(d / 'other_counts.csv.gz')
                z = execute_pair_model(d, comparator, ep, rscript)
                z['variant'] = variant; z['comparator'] = comparator
                pair_frames.append(z)
            aa = csv(parent / 'unfiltered_matched/background.csv')
            bb = csv(parent / 'filtered_matched/background.csv')
            require(aa.gene_symbol.equals(bb.gene_symbol) and np.array_equal(boolean(aa.retained, 'background'), boolean(bb.retained, 'background')), 'Background differs between matched variants')
            frames.extend(pair_frames)
            statuses.append(dict(comparator=comparator, status='completed_matched_sensitivity', donors=pairs.donor.nunique(),
                                 common_genes=int(aa.retained.sum()), reason='shared labels;fixed automatic filter;not doublet truth'))
        except NotTestable as e:
            statuses.append(dict(comparator=comparator, status='not_testable', reason=str(e)))
        except Exception as e:
            statuses.append(dict(comparator=comparator, status='technical_error', reason=str(e)))
            (parent / 'ERROR.txt').write_text(traceback.format_exc(), encoding='utf-8')
        save(out / '298_model_status.csv', statuses)
    grid = fixed_grid(dict(endpoints=ep, frames=frames))
    status_table = pd.DataFrame(statuses)[['comparator', 'status', 'reason']]
    grid = grid.merge(status_table, on='comparator', how='left', validate='many_to_one')
    save(out / '298_all_fixed_interaction_results.csv', grid)
    summary = grid.groupby(['variant', 'comparator'], sort=False).agg(planned=('gene_symbol', 'size'),
        tested=('testable', 'sum'), q1422_below005=('q_1422_within_variant', lambda s: int(s.lt(.05).sum()))).reset_index()
    save(out / '298_interaction_summary.csv', summary)
    un = grid[grid.variant == 'unfiltered_matched'].drop(columns='variant')
    fi = grid[grid.variant == 'filtered_matched'].drop(columns='variant')
    paired = un.merge(fi, on=['comparator', 'gene_symbol'], suffixes=('_unfiltered', '_filtered'), validate='one_to_one')
    paired['both_testable'] = paired.testable_unfiltered & paired.testable_filtered
    paired['coefficient_delta_descriptive'] = paired.logFC_filtered - paired.logFC_unfiltered
    paired['direction_same'] = np.where(paired.both_testable, np.sign(paired.logFC_unfiltered) == np.sign(paired.logFC_filtered), None)
    paired['difference_hypothesis_test_performed'] = False
    save(out / '299_all_paired_interactions.csv', paired)
    save(out / '299_focal_gene_results.csv', grid[grid.gene_symbol.isin(FOCAL)])
    save(out / '299_fixed_endpoints.csv', ep)
    return statuses
