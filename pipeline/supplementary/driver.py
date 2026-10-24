"""Supplementary analyses from this run's public inputs, without result imports."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gc
import hashlib
import importlib
import importlib.metadata
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd

from pipeline.io import dump, require, sha
from pipeline.lifecycle import claim, code_identity, assert_code, current_source, parent_run
from pipeline.matched import lock_output
from .common import ROOT, RDIR, Download, CurrentFiles, readj, save, reserved_bh, discovery_direction_comparison

GROUPS = ('age', 'quality', 'bulk', 'cognition', 'external', 'astro', 'interactions', 'mixtures', 'simulations')
STOP = threading.Event()


def check_python(groups):
    from packaging.version import Version
    modules = {'numpy': 'numpy', 'pandas': 'pandas', 'scipy': 'scipy', 'h5py': 'h5py'}
    if 'external' in groups:
        modules['openpyxl'] = 'openpyxl'
    if set(groups) & {'astro', 'interactions'}:
        modules.update({'scikit-learn': 'sklearn', 'threadpoolctl': 'threadpoolctl'})
    if 'interactions' in groups:
        modules['scrublet'] = 'scrublet'
    versions, failures = {}, []
    for package, module in modules.items():
        try:
            importlib.import_module(module)
            versions[package] = importlib.metadata.version(package)
            bounds = {'scikit-learn': ('1.4', '2'), 'threadpoolctl': ('3', '4'),
                      'scrublet': ('0.2.3', '0.3')}
            if package in bounds:
                lo, hi = bounds[package]
                require(Version(lo) <= Version(versions[package]) < Version(hi),
                        f'{package} {versions[package]} is outside [{lo}, {hi})')
        except Exception as exc:
            failures.append(f'{package}: {exc}')
    require(not failures, 'Supplementary Python dependencies are not loadable:\n' + '\n'.join(failures) +
            '\nNo packages were downloaded. Select an existing compatible environment.')
    return versions


def rrun(args, script, arguments, log):
    require(not STOP.is_set(), 'Execution was interrupted')
    env = os.environ.copy()
    libraries = [str(p.resolve()) for p in args.r_library]
    if env.get('R_LIBS'):
        libraries.extend(env['R_LIBS'].split(os.pathsep))
    env['R_LIBS'] = os.pathsep.join(dict.fromkeys(libraries))
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[key] = '1'
    with Path(log).open('w', encoding='utf-8') as stream:
        proc = subprocess.run([args.rscript, '--vanilla', str(RDIR / script), *map(str, arguments)],
                              stdout=stream, stderr=subprocess.STDOUT, env=env)
    require(proc.returncode == 0, f'R failed ({proc.returncode}); see {log}')


def panel(core):
    frame = pd.read_csv(core / 'results/panel.csv')
    require(len(frame) and frame.gene_symbol.is_unique and np.isfinite(frame.beta).all(), 'Invalid current panel')
    frame['discovery_direction'] = np.where(frame.beta > 0, 'Up', 'Down')
    require(frame.beta.ne(0).all(), 'Zero-effect discovery endpoint')
    return frame


def discovery_endpoints(core):
    full = pd.read_csv(core / 'results/panel_discovery.csv.gz')
    selected = full[full.q.lt(.05)].copy()
    selected['discovery_region'] = selected.region.str.replace(' ', '_', regex=False)
    selected['discovery_direction'] = np.where(selected.beta > 0, 'Up', 'Down')
    selected['discovery_log2_fold_change'] = selected.beta
    require(not selected.duplicated(['discovery_region', 'gene_symbol']).any(), 'Duplicate discovery slots')
    return selected


def age(core, out, dl, args):
    inp = out / 'inputs'; inp.mkdir()
    m = pd.read_csv(core / 'prepared/discovery_metadata.csv')
    m = m.rename(columns={'unit_id': 'analysis_unit_id', 'donor_id': 'parsed_individual_id',
        'region': 'parsed_brain_region', 'group': 'parsed_group', 'age': 'parsed_age', 'sex': 'parsed_sex'})
    m['parsed_brain_region'] = m.parsed_brain_region.str.replace(' ', '_', regex=False)
    m['parsed_group'] = m.parsed_group.replace({'ad': 'AD'})
    m.to_csv(inp / 'metadata.csv', index=False)
    g = pd.read_csv(core / 'results/panel_discovery.csv.gz').rename(columns={'beta': 'beta_input_scale', 'q': 'BH_q'})
    g['region'] = g.region.str.replace(' ', '_', regex=False); g['model'] = 'linear_age'
    g.to_csv(inp / 'models.csv', index=False)
    discovery_endpoints(core).to_csv(inp / 'fixed.csv', index=False)
    shutil.copyfile(core / 'prepared/discovery_expression.csv.gz', inp / 'expression.csv.gz')
    rrun(args, 'age.R', [inp, out], out / 'engine.log')


def bulk(core, out, dl, args):
    from . import bulk as B
    files = {k: dl.get(url, 'GSE48350_' + k + ('.txt.gz' if k == 'matrix' else '.soft'))
             for k, url in B.URLS.items()}
    rows, selected = B.parse_metadata(files['metadata'], out)
    B.annotation(files['platform'], out)
    B.extract_matrix(files['matrix'], rows, selected, out)
    discovery_endpoints(core).to_csv(out / '186_fixed_endpoints.csv', index=False)
    rrun(args, 'bulk.R', [out], out / 'engine.log')


def quality(core, out, dl, args):
    shutil.copyfile(core / 'prepared/human/metadata.csv', out / 'metadata.csv')
    discovery_endpoints(core).to_csv(out / 'endpoints.csv', index=False)
    rrun(args, 'quality.R', [str(out).encode().hex(), str(core / 'prepared/human').encode().hex()], out / 'engine.log')
    z = pd.read_csv(out / 'quality_variants.csv')
    require(not z.status.astype(str).str.contains('error', case=False).any(), 'Quality variant fitting error')
    z['significant'] = z.sensitivity_q.lt(.05)
    z.groupby(['variant', 'arm', 'model'], as_index=False).agg(
        planned=('gene_symbol', 'size'), significant=('significant', 'sum')
    ).to_csv(out / 'quality_model_variants.csv', index=False)


def cognition(core, out, dl, args):
    from . import cognition as C
    files = {name: dl.get(url, name) for name, url in C.SOURCES.items()}
    C.prepare(files['GSE1297_series_matrix.txt.gz'], files['GPL96.annot.gz'], out)
    panel(core).to_csv(out / '229_fixed_endpoints.csv', index=False)
    rrun(args, 'cognition.R', [str(out).encode('utf-8').hex()], out / 'engine.log')


def external(core, out, dl, args):
    from . import external as E
    fixed = panel(core).rename(columns={'gene_symbol': 'human_gene_symbol'})
    orth = pd.read_csv(core / 'results/ortholog_panel.csv', keep_default_na=False)
    fixed = fixed.merge(orth[['gene_symbol', 'mouse_gene']].rename(columns={'gene_symbol': 'human_gene_symbol'}),
                        on='human_gene_symbol', how='left', validate='one_to_one')
    fixed['gene_symbol'] = fixed.mouse_gene.fillna('')
    fixed['ortholog_mapped'] = fixed.gene_symbol.ne('')
    options = SimpleNamespace(rscript=args.rscript, mouse_manifest=core / 'prepared/mouse_units.csv',
                              human_cell_labels=None)
    rows = []
    for step, accession in E.STUDIES.items():
        d = out / accession; d.mkdir()
        z, state = getattr(E, 'run' + str(step))(d, dl, fixed, options)
        z = E.adjust_results(z)
        z['source_accession'] = accession; z['step'] = step
        z = discovery_direction_comparison(z)
        z.to_csv(d / 'all_fixed_results.csv', index=False)
        dump(d / 'status.json', {'status': state})
        rows.append(z)
    z = pd.concat(rows, ignore_index=True)
    ix = z.primary_task.eq(True) & z.human_gene_symbol.isin(E.FOCUS)
    require(int(ix.sum()) == 36 and not z.loc[ix].duplicated(['source_accession', 'task', 'human_gene_symbol']).any(),
            'Primary focal family must retain exactly 36 slots')
    z['primary_focus_36_member'] = ix
    z['primary_family_size'] = np.where(ix, 36, np.nan)
    z['q_across_primary_directions'] = np.nan
    z.loc[ix, 'q_across_primary_directions'] = reserved_bh(z.loc[ix, 'P_Value'])
    z.to_csv(out / 'all_fixed_endpoint_results.csv', index=False)
    z[z.human_gene_symbol.isin(E.FOCUS)].to_csv(out / 'focal_results.csv', index=False)


def astro(core, out, dl, args):
    from . import astro as A
    readme = dl.get(A.geo(A.REFERENCE) + 'suppl/GSE243292_processed_data_readme.txt')
    require('Count matrix' in readme.read_text(encoding='utf-8-sig'), 'Reference source count declaration changed')
    register = A.target_register(A.soft(dl.get(A.geo(A.TARGET) + 'soft/' + A.TARGET + '_family.soft.gz')))
    register.to_csv(out / 'target_register.csv', index=False)
    fixed = panel(core).gene_symbol.tolist()
    allowed = None
    for url in register['features.tsv.gz']:
        f = pd.read_csv(dl.get(url), sep='\t', header=None, dtype=str, keep_default_na=False)
        genes = set(f.loc[f[2].eq('Gene Expression'), 1] if f.shape[1] >= 3 else f[1])
        allowed = genes if allowed is None else allowed & genes
    x, genes, obs = A.get_reference(dl)
    model, features, ast, oli, gates = A.classifier(x, genes, obs, set(fixed), out, allowed=allowed)
    del x; gc.collect()
    counts, units, census = A.prepare_target(register, dl, model, features, out)
    first = out / 'external_models'; first.mkdir()
    second = out / 'composition_state'; second.mkdir()
    A.external_models(first, counts, units, ast, fixed, args.rscript)
    A.composition_and_state(second, counts, units, census, obs, ast, oli, gates, fixed, args.rscript)
    # Missing annotation eligibility is reported explicitly, never filled with p=0.
    focus = pd.concat([pd.read_csv(first / '272_focal_results.csv').assign(module='external'),
                       pd.read_csv(second / '273_focal_state_context_results.csv').assign(module='state')], ignore_index=True)
    primary = focus.task.eq('astro_primary_all_donors')
    focus['reporting_family_size'] = np.where(primary, 3, 12)
    focus['reporting_BH'] = np.nan
    for mask, size in ((primary, 3), (~primary, 12)):
        values = focus.loc[mask, 'P_Value'].to_numpy(float)
        require(len(values) <= size, 'Unexpected focal-family size')
        focus.loc[mask, 'reporting_BH'] = reserved_bh(np.r_[values, np.full(size-len(values), np.nan)])[:len(values)]
    focus.to_csv(out / 'focal_results.csv', index=False)


def interactions(core, out, dl, args):
    from . import annotation as A, interactions as I, interaction_math as U
    samples = A.target_register(A.soft(dl.get(A.geo(A.TARGET) + 'soft/' + A.TARGET + '_family.soft.gz')))
    samples.to_csv(out / 'target_register.csv', index=False)
    paths = {row.gsm: dl.get(row.count_url) for row in samples.itertuples(index=False)}
    fixed = panel(core)[['gene_symbol', 'discovery_direction']]
    identity = A.obtain_identity(dl, out, out)
    original, units, ast, gates, mapped = A.annotate_and_aggregate(samples, paths, dl, fixed, out, identity)
    readiness = []
    for fraction, cell_type in itertools.product(('P', 'V'), ('Ast', 'End')):
        z = units[units.fraction.eq(fraction) & units.cell_type.eq(cell_type)].copy()
        z['cell_count_eligible'] = z.n_cells.ge(U.POLICY['min_cells'])
        complete = z[['age', 'sex', 'pmi', 'ph_cerebellum']].notna().all(axis=1)
        eligible = z[z.cell_count_eligible & complete]
        readiness.append(dict(fraction=fraction, cell_type=cell_type, units=len(z),
            min_cells=U.POLICY['min_cells'], eligible_donors=eligible.donor.nunique(),
            eligible_AD=int(eligible.case.eq(1).sum()), eligible_control=int(eligible.case.eq(0).sum()),
            missing_covariates=int((~complete).sum())))
    save(out / 'cell_class_readiness.csv', readiness)
    primary = units[(units.fraction == 'P') & units.cell_type.eq(ast) & units.n_cells.ge(50)].copy()
    primary = primary.dropna(subset=['age', 'sex', 'pmi', 'ph_cerebellum'])
    d = out / '288_P_models'; d.mkdir(); primary.to_csv(d / 'metadata.csv', index=False)
    vectors, metadata, risks, repeats = {}, [], [], []
    for sample in samples.itertuples(index=False):
        print('DOUBLET', sample.gsm, flush=True)
        x, genes, barcodes, mapping = A.read_counts(paths[sample.gsm])
        ann = pd.read_csv(out / f'285_{sample.gsm}_cell_annotations.csv.gz')
        require(np.array_equal(ann.barcode, barcodes), 'Fresh annotation barcode order changed')
        keep, _, _, _ = A.qc_counts(x, genes)
        require(np.array_equal(I.boolean(ann.qc_pass, 'qc_pass'), keep), 'Annotation and doublet QC disagree')
        d = out / 'doublet_libraries' / sample.gsm
        retained, row = U.doublet_risk(x, barcodes, keep, d, sample, ann)
        risks.append(row)
        save(out / 'doublet_library_summary.csv', risks)
        if retained is None:
            require(row['status'] == 'not_testable_insufficient_qc_cells' and int(keep.sum()) < 100,
                    'Unexpected missing doublet decision')
            print('DOUBLET NOT TESTABLE', sample.gsm, 'QC cells:', int(keep.sum()), flush=True)
            continue
        if row['threshold_review_required']:
            observed = pd.read_csv(d / 'observed_scores.csv.gz')
            repeats.extend(I.seed_diagnostic(x, keep, observed, sample, out / 'doublet_seeds' / sample.gsm))
            require(not any(r['status'] == 'technical_error' for r in repeats), 'Doublet seed diagnostic failed')
        for unit in units[units.gsm.eq(sample.gsm)].itertuples(index=False):
            mask = retained & ann.label.eq(unit.cell_type).to_numpy()
            values = np.asarray(x[mask].sum(axis=0)).ravel().astype(np.int64)
            require(np.array_equal(genes, original.index), 'Filtered gene universe changed')
            vectors[unit.unit_id] = values
            rec = unit._asdict(); rec.update(n_cells=int(mask.sum()), total_umi=int(values.sum()))
            metadata.append(rec)
        save(out / 'doublet_seed_summary.csv', repeats, columns=None if repeats else
             ['gsm', 'seed', 'status', 'threshold', 'n_retained'])
    save(out / 'doublet_seed_summary.csv', repeats, columns=None if repeats else
         ['gsm', 'seed', 'status', 'threshold', 'n_retained'])
    filtered = pd.DataFrame(vectors, index=original.index)
    filtered.index.name = 'gene_symbol'; filtered.to_csv(out / 'filtered_counts.csv.gz')
    units = save(out / 'filtered_units.csv', metadata, columns=list(units.columns))
    statuses = I.models(CurrentFiles(out), out, original, filtered, units, args.rscript)
    require(not any(row['status'] == 'technical_error' for row in statuses), 'Interaction model execution failed')


def mixtures(core, out, dl, args):
    aux = out / 'inputs'; aux.mkdir()
    shutil.copyfile(core / 'prepared/human/metadata.csv', aux / 'sea_metadata.csv')
    panel(core).to_csv(aux / 'panel.csv', index=False)
    rrun(args, 'mixtures.R', ['cell', core / 'prepared/human/expression.csv.gz', out, aux], out / 'engine.log')


def animal_metadata(core):
    m = pd.read_csv(core / 'prepared/mouse_units.csv', dtype={'cage': str, 'slide': str})
    require(m.genotype.isin(['WT', '5X']).all() and m.sex.isin(['M', 'F']).all(), 'Unknown mouse genotype/sex coding')
    m['genotype'] = m.genotype.eq('5X').astype(int); m['male'] = m.sex.eq('M').astype(int)
    return m


def record_simulation_order(core, out, actual):
    """Freeze random-draw positions without changing the scientific input order."""
    fields = ['SAMPLE', 'cage', 'genotype', 'age', 'male', 'slide']
    require(set(fields) <= set(actual), 'Simulation animal identifiers/covariates missing')
    require(not actual[fields].isna().any().any() and not actual.SAMPLE.duplicated().any(),
            'Missing or duplicate simulation animal identifiers')
    ordered = actual[fields].copy().reset_index(drop=True)
    ordered.insert(0, 'draw_position', np.arange(1, len(ordered) + 1))
    content = ordered.to_csv(index=False, lineterminator='\n').encode('utf-8')
    source = core / 'prepared/mouse_units.csv'
    record = dict(schema='simulation_animal_order_v1', rows=len(ordered),
                  input='core/prepared/mouse_units.csv', input_sha256=sha(source),
                  order_sha256=hashlib.sha256(content).hexdigest(),
                  order_policy='preserve_prepared_input_order',
                  two_model_seed=319324, three_model_seed_base=3252026)
    path = out / 'simulation_input_order.json'
    csv_path = out / 'simulation_animal_order.csv'
    if path.exists():
        require(readj(path) == record, 'Simulation input/order changed; use a fresh experiment')
        require(csv_path.is_file() and sha(csv_path) == record['order_sha256'],
                'Simulation animal-order record is missing or damaged')
    else:
        require(not csv_path.exists() or sha(csv_path) == record['order_sha256'],
                'Existing simulation animal-order record differs')
        csv_path.write_bytes(content)
        dump(path, record)
    return record


def simulations(core, out, dl, args):
    from . import layouts as L, conditional_simulation as C
    actual = animal_metadata(core)
    record_simulation_order(core, out, actual)
    layouts, registry = L.build_layouts(3252026, actual)
    registry.to_csv(out / 'layout_registry.csv', index=False)
    pd.concat([u.assign(layout_id=name) for name, (u, _, _) in layouts.items()],
              ignore_index=True).to_csv(out / 'layout_units.csv', index=False)
    designs = {}
    for scope in ('whole', 'cortex'):
        for kind in ('linear', 'age_factor', 'age_factor_slide'):
            u, b, f = layouts['actual_' + kind]
            # This simulation estimates the second coefficient in both designs.
            order = [c for c in f if c != 'genotype']; order.insert(1, 'genotype')
            designs[scope + ':' + kind] = (b.to_numpy(float), f[order].to_numpy(float), u.cage.to_numpy())
    C.simulate_designs(designs, out, 1737, 319324)
    rrun(args, 'layouts.R', ['preflight', out / 'R_versions.csv', out / 'R_session.txt'], out / 'preflight.log')
    scenarios = []
    for name, (u, b, f) in layouts.items():
        d = out / name; d.mkdir(exist_ok=True)
        u.to_csv(d / 'units.csv', index=False); b.to_csv(d / 'base.csv', index=False); f.to_csv(d / 'fixed.csv', index=False)
        for rho, shift, truth in itertools.product((0., .05, .20), (0., .5), (0., .3)):
            sid = f'{name}_rho{rho:g}_shift{shift:g}_beta{truth:g}'
            seed = int(hashlib.sha256(f'3252026:{sid}'.encode()).hexdigest()[:8], 16) % 1_000_000_000
            rec = dict(scenario_id=sid, layout_id=name, rho=rho, shift=shift, truth=truth, seed=seed,
                units=str(d / 'units.csv'), base_design=str(d / 'base.csv'), fixed_design=str(d / 'fixed.csv'),
                output=str(d / (sid + '.csv.gz')), background_genes=64, reps=625,
                fixed_estimable=bool(registry.set_index('layout_id').loc[name, 'fixed_estimable']))
            scenarios.append(rec)
    schedule = pd.DataFrame(scenarios); schedule.to_csv(out / 'scenarios.csv', index=False)
    def execute(rec):
        target = Path(rec['output']); receipt = target.with_suffix('.complete.json')
        identity = hashlib.sha256(json.dumps(rec, sort_keys=True).encode()).hexdigest()
        if receipt.exists():
            old = readj(receipt)
            require(old['identity'] == identity and target.is_file() and old['sha256'] == sha(target), 'Damaged scenario result')
        else:
            config = target.with_suffix('.config.csv'); save(config, [rec])
            rrun(args, 'layouts.R', ['simulate', config], target.with_suffix('.log'))
            z = pd.read_csv(target)
            require(len(z) == 3 * rec['reps'] and z.scenario_id.eq(rec['scenario_id']).all()
                    and not z.duplicated(['replicate', 'model']).any(), 'Incomplete layout replicates')
            dump(receipt, {'identity': identity, 'sha256': sha(target)})
        return target
    frames = []
    pool = ThreadPoolExecutor(max_workers=args.workers)
    try:
        futures = [pool.submit(execute, rec) for rec in scenarios]
        for i, future in enumerate(as_completed(futures), 1):
            frames.append(pd.read_csv(future.result()))
            print(f'LAYOUTS {i}/{len(futures)}', flush=True)
    except BaseException:
        STOP.set()
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    else:
        pool.shutdown(wait=True)
    raw = pd.concat(frames, ignore_index=True)
    raw.to_csv(out / 'all_replicates.csv.gz', index=False)
    summary, paired = L.summarize(raw, schedule[['scenario_id','layout_id','rho','shift','truth']], registry)
    summary.to_csv(out / 'model_summary.csv', index=False)
    paired.to_csv(out / 'paired_model_comparisons.csv', index=False)
    raw[raw.status.ne('ok')].to_csv(out / 'unavailable_fits.csv.gz', index=False)


def run_stage(name, root, callback):
    receipt = root / (name + '.complete.json')
    if receipt.exists():
        info = readj(receipt)
        directory = Path(info['directory']).resolve()
        require(directory.is_relative_to(root / name) and directory.is_dir(), 'Invalid current-run stage directory')
        for relative, expected in info['files'].items():
            path = (root / relative).resolve()
            require(path.is_relative_to(directory) and path.is_file() and sha(path) == expected,
                    'Incomplete/damaged same-run stage: ' + relative)
        print('RESUME CURRENT RUN', name, flush=True)
        return directory
    attempts = root / name; attempts.mkdir(exist_ok=True)
    # Keep interrupted files separate; never adopt a partial fit as completed work.
    if name == 'simulations' and (attempts / 'attempt_001').exists():
        out = attempts / 'attempt_001'
    else:
        number = len(list(attempts.glob('attempt_*'))) + 1
        out = attempts / f'attempt_{number:03d}'; out.mkdir()
    callback(out)
    files = {p.relative_to(root).as_posix(): sha(p) for p in out.rglob('*') if p.is_file()}
    require(files, 'Stage produced no output: ' + name)
    dump(receipt, {'directory': str(out), 'files': files, 'completed': time.time()})
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--core', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--groups', nargs='+', choices=GROUPS, default=list(GROUPS))
    p.add_argument('--raw-cache', type=Path, action='append', default=[])
    p.add_argument('--r-library', type=Path, action='append', default=[])
    p.add_argument('--rscript', default=shutil.which('Rscript') or 'Rscript')
    p.add_argument('--workers', type=int, default=60)
    p.add_argument('--resume', action='store_true')
    p.add_argument('--offline', action='store_true')
    args = p.parse_args(argv)
    require(args.workers > 0, 'Positive worker count required')
    versions = check_python(args.groups)
    core = args.core.resolve(); root = args.output.resolve()
    run_id = current_source(core, parent_run())
    require(core.parent == root.parent, 'Supplement must be adjacent to this run core')
    root.mkdir(parents=True, exist_ok=True)
    identity = code_identity(ROOT)
    with lock_output(root):
        claim(root, 'supplementary', args.resume, run_id, identity)
        software = root / 'python.json'
        if software.exists():
            require(readj(software) == versions, 'Supplementary package versions changed during this run')
        else:
            dump(software, versions)
        library = [str(x.resolve()) for x in args.r_library] + os.environ.get('R_LIBS', '').split(os.pathsep)
        os.environ['R_LIBS'] = os.pathsep.join(dict.fromkeys(x for x in library if x))
        dl = Download(root / 'raw', args.raw_cache, args.offline)
        state = {'status': 'RUNNING', 'groups': args.groups, 'run_id': run_id, 'workers': args.workers,
                 'started': time.time(), 'completed_groups': {}, 'pid': os.getpid()}
        dump(root / 'status.json', state)
        try:
            for name in args.groups:
                assert_code(ROOT, identity)
                require(check_python(args.groups) == versions, 'Supplementary package versions changed during execution')
                state['current_group'] = name; dump(root / 'status.json', state)
                result = run_stage(name, root, lambda out: globals()[name](core, out, dl, args))
                state['completed_groups'][name] = str(result)
                dump(root / 'status.json', state)
            assert_code(ROOT, identity)
            require(check_python(args.groups) == versions, 'Supplementary package versions changed during execution')
            state['status'] = 'REQUESTED_SUPPLEMENTARY_GROUPS_COMPLETE'
        except BaseException as exc:
            state.update(status='INTERRUPTED' if isinstance(exc, KeyboardInterrupt) else 'FAILED', reason=str(exc))
            raise
        finally:
            state['updated'] = time.time(); dump(root / 'status.json', state)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
