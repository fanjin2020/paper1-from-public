#!/usr/bin/env python3
"""Collect completed Paper1 evidence into a verified directory, never a ZIP.

Standard library only. Does not run analyses, install packages, modify source,
or write inside the experiment. Run on the server after reporting completes.
"""
import argparse
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import sys
import time


def default_review_output(root, now=None):
    """Keep transferred paths short, independent of the experiment name."""
    return root.with_name('p1_review_' + (now or datetime.now()).strftime('%Y%m%d_%H%M%S'))


DEFAULT_PROJECT = Path.cwd()
DEFAULT_RUN = DEFAULT_PROJECT.parent / 'paper1_complete_from_public_20261022_v2'
GROUPS = ('age', 'quality', 'bulk', 'cognition', 'external', 'astro',
          'interactions', 'mixtures', 'simulations')
STAGES = ('download', 'prepare_panel', 'panel', 'orthologs', 'prepare_human',
          'prepare_mouse', 'human', 'mouse', 'calibration')
TABLE_SUFFIXES = ('.csv', '.csv.gz', '.tsv', '.tsv.gz')
FINAL_SUFFIXES = TABLE_SUFFIXES + ('.json', '.txt', '.md', '.pdf', '.png', '.xlsx')


def need(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    with path.open(encoding='utf-8-sig') as handle:
        return json.load(handle)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=True) + '\n', encoding='utf-8')


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024**2), b''):
            digest.update(chunk)
    return digest.hexdigest()


def inside(root, relative):
    """Accept portable relative paths, rejecting symlinks and parent traversal."""
    value = str(relative).replace('\\', '/')
    p = PurePosixPath(value)
    need(value and not p.is_absolute() and '..' not in p.parts
         and ':' not in value, 'Unsafe relative path: ' + value)
    candidate = root.joinpath(*p.parts)
    for parent in (candidate, *candidate.parents):
        if parent == root:
            break
        need(not parent.is_symlink(), 'Symlink not collected: ' + str(parent))
    need(candidate.resolve().is_relative_to(root.resolve()), 'Path escapes root: ' + value)
    return candidate


def csv_records(path):
    opener = gzip.open if path.name.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        yield from csv.DictReader(handle)


def table_info(path):
    if not path.name.endswith(TABLE_SUFFIXES):
        return '', ''
    opener = gzip.open if path.name.endswith('.gz') else open
    delimiter = '\t' if '.tsv' in path.name else ','
    previous=csv.field_size_limit()
    csv.field_size_limit(max(previous,2**31-1))
    try:
        with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
            rows = csv.reader(handle, delimiter=delimiter)
            header = next(rows, [])
            return sum(1 for _ in rows), len(header)
    finally:
        csv.field_size_limit(previous)


def write_csv(path, records, fields):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


@contextmanager
def read_locks(root):
    handles = []
    try:
        if sys.platform != 'win32':
            import fcntl
            for folder in ('', 'core', 'matched', 'benchmark', 'supplementary', 'reporting'):
                path = root / folder / 'run.lock'
                if path.exists():
                    handle = path.open('rb')
                    handles.append(handle)
                    try:
                        fcntl.flock(handle.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
                    except BlockingIOError as exc:
                        raise ValueError('Writer still active; collect after completion: ' + str(path)) from exc
        yield
    finally:
        for handle in reversed(handles):
            handle.close()


class Plan:
    def __init__(self, root):
        self.root = root
        self.items = {}
        self.expected = {}
        self.warnings = []
        self.excluded = []
        self.checks = {}

    def expect(self, source, digest, size=None):
        need(len(digest) == 64, 'Invalid receipt digest: ' + str(source))
        key = str(source)
        old = self.expected.get(key)
        need(not old or old[0] == digest, 'Conflicting completion hashes: ' + str(source))
        self.expected[key] = (digest, size if size is not None else (old[1] if old else None))

    def add(self, relative, category='completion', required=True):
        source = inside(self.root, relative)
        if not source.is_file():
            need(not required, 'Missing required file: ' + str(source))
            return
        self.items[source.relative_to(self.root).as_posix()] = (source, category)

    def tree(self, relative, category):
        folder = inside(self.root, relative)
        need(folder.is_dir(), 'Missing result directory: ' + str(folder))
        for source in sorted(folder.rglob('*')):
            need(not source.is_symlink(), 'Symlink not collected: ' + str(source))
            if source.is_file() and source.name.endswith(FINAL_SUFFIXES):
                self.add(source.relative_to(self.root), category)

    def source(self, source, relative):
        need(not source.is_symlink() and source.is_file(), 'Invalid source file: ' + str(source))
        self.items['source_snapshot/' + relative] = (source, 'source_snapshot')


def validate_completion(plan):
    root = plan.root
    state = read_json(root / 'state.json')
    report = read_json(root / 'reporting/status.json')
    need(state.get('status') == 'INTEGRATED_ANALYSES_COMPLETE',
         'Run not complete: ' + str(state.get('status')) +
         '. Complete the reporting stage before collecting; see docs/REEXPORT.md.')
    need(report.get('status') == 'RESULT_EXPORTS_COMPLETE', 'Final reporting has not completed')
    run_id = state.get('run_id')
    need(run_id and report.get('run_id') == run_id, 'Run/report identities differ')
    identities = []
    for folder in ('', 'core', 'matched', 'benchmark', 'supplementary'):
        relative = str(PurePosixPath(folder) / 'fresh_run.json')
        marker = read_json(root / relative)
        need(marker.get('run_id') == run_id and marker.get('policy') == 'raw_downloads_only_v1',
             'Mixed or unsupported run identity: ' + relative)
        identities.append(marker.get('code'))
        plan.add(relative)
    need(identities[0] and len(set(identities)) == 1, 'Component source identities differ')
    core = read_json(root / 'core/status.json')
    matched = read_json(root / 'matched/state.json')
    summary = read_json(root / 'matched/results/analysis_summary.json')
    benchmark = read_json(root / 'benchmark/run_state.json')
    supp = read_json(root / 'supplementary/status.json')
    need(core.get('status') == 'completed', 'Core analysis incomplete')
    need(matched.get('status') == 'COMPLETED' and not matched.get('validation_only'),
         'Matched analysis incomplete or smoke-only')
    need(summary.get('status') == 'ANALYSES_COMPLETE' and not summary.get('validation_only'),
         'Matched results require review or are smoke-only')
    need(benchmark.get('status') in ('COMPLETED', 'COMPLETED_WITH_FIT_FAILURES'),
         'Benchmark incomplete')
    need(supp.get('status') == 'REQUESTED_SUPPLEMENTARY_GROUPS_COMPLETE'
         and supp.get('run_id') == run_id, 'Supplementary analyses incomplete')
    need(set(supp.get('completed_groups', {})) == set(GROUPS), 'Not all nine supplementary groups completed')
    if benchmark['status'] == 'COMPLETED_WITH_FIT_FAILURES':
        plan.warnings.append('Benchmark records fit failures; inspect availability/failure tables.')
    plan.checks.update(run_id=run_id, source_identity=identities[0],
                       manuscript_agreement='NOT_ASSESSED')


def add_reporting(plan):
    root = plan.root
    records = list(csv_records(root / 'reporting/source_index.csv'))
    need({r['section'] for r in records} == {f'S{i:02}' for i in range(1, 21)},
         'Expected all 20 supplementary sections')
    files = {r['file'] for r in records}
    need(len(files) == len(records), 'Duplicate source-index paths')
    for record in records:
        source = inside(root / 'reporting', record['file'])
        plan.expect(source, record['sha256'], int(record['bytes']))
        plan.add(source.relative_to(root), 'supplementary_table')
    mapping = list(csv_records(root / 'reporting/block_sources.csv'))
    declared = {(b, r['file']) for r in records for b in r['blocks'].split(',') if b}
    need({(r['block'], r['file']) for r in mapping} == declared, 'Table-block map disagrees with source index')
    need(len({r['block'] for r in mapping}) == 105, 'Expected 105 mapped table blocks')
    plan.tree('reporting', 'reporting')
    primary = list(csv_records(root / 'figures/figure_index.csv'))
    secondary = list(csv_records(root / 'reporting/figures/figure_index.csv'))
    need(len(primary) == len({r['figure'] for r in primary}) == 8, 'Expected 8 main figures')
    need(len(secondary) == len({r['figure'] for r in secondary}) == 12, 'Expected 12 supplementary figures')
    for row in primary:
        for column in ('pdf', 'csv'):
            plan.add('figures/' + row[column], 'main_figure')
        plan.add('figures/' + row['figure'] + '.png', 'main_figure')
    for row in secondary:
        need(row['numeric_source'] in files, 'Supplementary figure source not indexed: ' + row['figure'])
        for suffix in ('.pdf', '.png', '.csv'):
            plan.add('reporting/figures/' + row['figure'] + suffix, 'supplementary_figure')
    plan.tree('figures', 'main_figure')
    plan.checks.update(sections=20, table_blocks=105, main_figures=8,
                       supplementary_figures=12, indexed_tables=len(records))


def add_final_results(plan):
    root = plan.root
    owner = read_json(root / 'core/owner.json')
    for stage in STAGES:
        relative = 'core/stages/' + stage + '.json'
        info = read_json(root / relative)
        need(info.get('signature') == owner.get('signature') and info.get('outputs'),
             'Core stage identity/outputs missing: ' + stage)
        for name, digest in info['outputs'].items():
            plan.expect(inside(root / 'core', name), digest)
        plan.add(relative)
    plan.tree('core/results', 'core_results')
    for name in ('core/prepared/human/metadata.csv', 'core/prepared/mouse_units.csv'):
        plan.add(name, 'sample_metadata')
    plan.tree('core/prepared/designs', 'design_matrix')
    receipt = read_json(root / 'matched/summary_receipt.json')
    need(receipt.get('protocol_hash') == read_json(root / 'matched/state.json').get('protocol_hash'),
         'Matched summary/protocol identities differ')
    need(receipt.get('files'), 'Missing matched-result hashes')
    for name, digest in receipt['files'].items():
        source = inside(root / 'matched', name)
        need(name.replace('\\', '/').split('/')[0] in ('results', 'figures'),
             'Unexpected matched final-result path: ' + name)
        plan.expect(source, digest)
        plan.add(source.relative_to(root), 'matched_final')
    benchmark = root / 'benchmark/final_results'
    records = list(csv_records(benchmark / 'file_manifest.csv'))
    need(records, 'Empty benchmark manifest')
    for row in records:
        source = inside(benchmark, row['path'])
        plan.expect(source, row['sha256'], int(row['bytes']))
        plan.add(source.relative_to(root), 'benchmark_final')
    plan.add('benchmark/final_results/file_manifest.csv')
    for required in ('human_sensitivity_fits.csv.gz', 'human_sensitivity_paired_changes.csv.gz',
                     'mouse_seven_arm_background_fits.csv.gz', 'mouse_seven_arm_panel_fits.csv.gz',
                     'count_pipeline_per_dataset.csv', 'count_pipeline_paired_per_dataset.csv',
                     'count_pipeline_unavailable_genes.csv', 'count_pipeline_diagnostic_inventory.csv'):
        key = 'matched/results/' + required
        need(key in plan.items, 'Essential detailed result absent from receipt: ' + key)


def add_supplementary_details(plan):
    root = plan.root
    for group in GROUPS:
        name = 'supplementary/' + group + '.complete.json'
        info = read_json(root / name)
        directory = Path(info['directory']).resolve()
        need(directory.is_relative_to(root / 'supplementary' / group) and directory.is_dir(),
             'Invalid completed attempt directory: ' + str(directory))
        need(info.get('files'), 'Empty supplementary completion receipt: ' + group)
        plan.add(name)
        for relative, digest in info['files'].items():
            source = inside(root / 'supplementary', relative)
            need(source.is_relative_to(directory), 'Receipt crosses supplementary attempts')
            plan.expect(source, digest)
            # Only top-level result tables, not library-level cells or matrices.
            low = source.name.lower()
            wanted = (source.parent == directory and source.name.endswith(FINAL_SUFFIXES)
                      and not any(term in low for term in ('counts', 'expression', 'barcodes', 'cell_metadata')))
            if wanted:
                plan.add(source.relative_to(root), 'supplementary_detail')
        plan.excluded.append({'scope': str(directory.relative_to(root)),
                              'not_collected': 'Counts, cell-level files and nested fit/batch directories; final indexed tables retained.'})


def add_records(plan):
    root = plan.root
    required = ('state.json', 'software.json', 'paper_coverage.json', 'core/status.json',
                'core/owner.json', 'core/run_config.json', 'core/software.json',
                'matched/state.json', 'matched/protocol.json', 'matched/input_manifest.json',
                'matched/summary_receipt.json', 'benchmark/run_state.json', 'benchmark/protocol.json',
                'supplementary/python.json', 'supplementary/status.json')
    for name in required:
        plan.add(name)
    for folder in ('matched/runtime', 'benchmark/runtime'):
        path = root / folder
        if path.is_dir():
            for source in path.glob('*.json'):
                plan.add(source.relative_to(root), 'environment')
    for base in ('core/raw', 'supplementary/raw'):
        for source in sorted((root / base).rglob('*.download.json')):
            plan.add(source.relative_to(root), 'raw_download_receipt')
    plan.add('core/raw_cache_import.json', 'raw_download_receipt', required=False)
    for source in sorted((root / 'maintenance').rglob('*.json')):
        plan.add(source.relative_to(root), 'maintenance_record')
    plan.excluded.extend([
        {'scope': 'raw inputs / prepared expression matrices', 'not_collected': 'FASTQ, CEL, MTX, H5/H5AD, RDS, NPZ and raw count matrices.'},
        {'scope': 'simulation execution', 'not_collected': 'Per-batch fit objects, runtime libraries, caches, logs and unfinished attempts. Final per-dataset count metrics and saved truth examples ARE retained.'},
        {'scope': 'benchmark individual batches', 'not_collected': 'Batch files are omitted. Final benchmark summaries/geometry/failure records are retained; this package does not independently refit models.'},
    ])


def add_source(plan, project):
    """Mirror the runner's code-identity contract without importing its code."""
    paths = [project / n for n in ('run.py', 'config.json', 'requirements.txt')]
    for name in ('reproduce_from_public.py', 'requirements-supplementary.txt'):
        if (project / name).is_file():
            paths.append(project / name)
    for folder in ('pipeline', 'analysis', 'tests', 'inputs'):
        paths.extend(p for p in (project / folder).rglob('*') if p.suffix in ('.py', '.R', '.json', '.csv'))
    digest = hashlib.sha256()
    for path in sorted(paths):
        path = inside(project, path.relative_to(project))
        need(path.is_file(), 'Source snapshot file missing: ' + str(path))
        relative = path.relative_to(project).as_posix()
        data = path.read_bytes()
        digest.update(relative.encode())
        digest.update(data)
        plan.source(path, relative)
        plan.expect(path, hashlib.sha256(data).hexdigest(), len(data))
    need(digest.hexdigest() == plan.checks['source_identity'],
         'Source differs from this completed run. Do not replace source just to collect results. '
         'Point --project-root to the matching source, or use --no-source and retain this limitation.')
    for name in ('README.md', 'README_CN.md', 'CITATION.cff', 'environment.yml', 'run_full_reproduction.sh'):
        path = inside(project, name)
        if path.is_file():
            plan.source(path, name)
    for path in sorted((project / 'docs').rglob('*.md')):
        safe = inside(project, path.relative_to(project))
        plan.source(safe, path.relative_to(project).as_posix())
    plan.checks['source_snapshot'] = 'EXECUTION_SOURCE_IDENTITY_MATCHED'


def build_plan(root, project=None):
    plan = Plan(root)
    validate_completion(plan)
    add_reporting(plan)
    add_final_results(plan)
    add_supplementary_details(plan)
    add_records(plan)
    if project is not None:
        add_source(plan, project)
    else:
        plan.checks['source_snapshot'] = 'NOT_INCLUDED'
        plan.warnings.append('Execution source not included (--no-source).')
    # Take a pre-copy fingerprint of unreceipted metadata and figure files too.
    print('CHECK selected final files and completion hashes', flush=True)
    last = time.monotonic()
    for i, (name, (source, category)) in enumerate(sorted(plan.items.items()), 1):
        stat = source.stat()
        digest = sha256(source)
        expectation = plan.expected.get(str(source))
        if expectation:
            need(digest == expectation[0] and (expectation[1] is None or stat.st_size == expectation[1]),
                 'File differs from completion receipt: ' + str(source))
        plan.expected[str(source)] = (digest, stat.st_size)
        if time.monotonic() - last > 5:
            print(f'CHECK {i}/{len(plan.items)} {name}', flush=True)
            last = time.monotonic()
    return plan


def copy_bundle(plan, out):
    need(not out.exists(), 'Output already exists; use a new output directory: ' + str(out))
    out.mkdir(parents=True)
    progress = out / 'COLLECTION_IN_PROGRESS.json'
    write_json(progress, {'source_run': str(plan.root), 'status': 'COPYING'})
    manifest = []
    try:
        last = time.monotonic()
        for i, (relative, (source, category)) in enumerate(sorted(plan.items.items()), 1):
            target = inside(out, relative)
            target.parent.mkdir(parents=True, exist_ok=True)
            expected, size = plan.expected[str(source)]
            shutil.copyfile(source, target)
            need(target.stat().st_size == size and sha256(target) == expected,
                 'Copy changed or was incomplete: ' + relative)
            rows, columns = table_info(target)
            manifest.append({'file': relative, 'category': category, 'bytes': size,
                             'sha256': expected, 'rows': rows, 'columns': columns})
            if time.monotonic() - last > 5:
                print(f'COPY {i}/{len(plan.items)} {relative}', flush=True)
                last = time.monotonic()
        # Verify the stored row counts, using the original CSV precision unchanged.
        counts = {row['file']: row['rows'] for row in manifest}
        for row in csv_records(out / 'reporting/source_index.csv'):
            need(counts['reporting/' + row['file']] == int(row['rows']),
                 'Export row count differs from source index: ' + row['file'])
        notes = {
            'source_run': str(plan.root), 'checks': plan.checks, 'warnings': plan.warnings,
            'excluded': plan.excluded,
            'verification_scope': 'Completion identities and selected final-file hashes. Raw input contents and omitted simulation batches were not rehashed or refit.',
        }
        write_json(out / 'collection_details.json', notes)
        readme = (
            '# Paper1 verification materials\n\n'
            'This is a read-only result collection, not a new analysis or a claim of manuscript agreement.\n\n'
            '* reporting/S01-S20: all indexed numeric tables and the 105-block source map.\n'
            '* figures: eight main figures and numeric source data.\n'
            '* reporting/figures: twelve supplementary figures and numeric source data.\n'
            '* core/results: human/mouse gene-level fits and calibration outputs.\n'
            '* matched/results: human sensitivity fits, seven mouse arms, per-dataset count metrics, unavailable-gene details and saved truth examples.\n'
            '* benchmark/final_results: benchmark summaries, design geometry and failure records.\n'
            '* supplementary: completed-attempt result tables, doublet/library eligibility, missingness and completion records.\n'
            '* source_snapshot: execution source when identity-matched (see collection_details.json).\n'
            '* completion/environment/download receipts: identify the run and its input provenance, without transferring raw inputs.\n\n'
            'Existing CSV/CSV.GZ/PDF/PNG files are byte-for-byte copies. No rounding, row filtering, XLSX conversion or ZIP archive is performed.\n'
            'CSV.GZ files remain in their original individually compressed form.\n'
            'Numbers can be compared to the manuscript and its Excel tables; Excel layout itself requires the local submission files.\n'
            'Raw counts, cells, runtime environments and large per-batch fit objects are intentionally absent.\n'
            'Hashes establish file consistency, not statistical validity or equivalence to another run.\n\n'
            'Only transfer this directory after COLLECTION_COMPLETE.json appears.\n'
            'After transfer: python paper1_collect_verification.py --verify PATH_TO_THIS_DIRECTORY\n'
        )
        (out / 'README_REVIEW.md').write_text(readme, encoding='utf-8')
        helper = Path(__file__).resolve()
        shutil.copyfile(helper, out / 'paper1_collect_verification.py')
        for name in ('collection_details.json', 'README_REVIEW.md', 'paper1_collect_verification.py'):
            target = out / name
            manifest.append({'file': name, 'category': 'collection_record', 'bytes': target.stat().st_size,
                             'sha256': sha256(target), 'rows': '', 'columns': ''})
        write_csv(out / 'file_manifest.csv', manifest,
                  ['file', 'category', 'bytes', 'sha256', 'rows', 'columns'])
        completion = {'status': 'VERIFICATION_PACKAGE_COMPLETE', 'run_id': plan.checks['run_id'],
                      'created_utc': datetime.now(timezone.utc).isoformat(),
                      'files': len(manifest), 'payload_bytes': sum(r['bytes'] for r in manifest),
                      'manifest_sha256': sha256(out / 'file_manifest.csv'),
                      'manuscript_agreement': 'NOT_ASSESSED'}
        progress.unlink()
        write_json(out / 'COLLECTION_COMPLETE.json', completion)
        return completion
    except BaseException as exc:
        write_json(out / 'COLLECTION_FAILED.json', {'status': 'INCOMPLETE', 'reason': str(exc)})
        raise


def verify_bundle(root):
    need(not (root / 'COLLECTION_IN_PROGRESS.json').exists() and not (root / 'COLLECTION_FAILED.json').exists(),
         'Collection is partial or failed')
    info = read_json(root / 'COLLECTION_COMPLETE.json')
    need(info.get('status') == 'VERIFICATION_PACKAGE_COMPLETE', 'Missing successful collection status')
    need(sha256(root / 'file_manifest.csv') == info['manifest_sha256'], 'Collection manifest changed')
    records = list(csv_records(root / 'file_manifest.csv'))
    need(len(records) == len({r['file'] for r in records}) == info['files'], 'Manifest count/uniqueness mismatch')
    for i, row in enumerate(records, 1):
        path = inside(root, row['file'])
        need(path.is_file() and path.stat().st_size == int(row['bytes']) and sha256(path) == row['sha256'],
             'Missing/changed transferred file: ' + row['file'])
        if i % 250 == 0:
            print(f'VERIFY {i}/{len(records)}', flush=True)
    expected = {r['file'] for r in records} | {'file_manifest.csv', 'COLLECTION_COMPLETE.json'}
    actual = {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
    need(actual == expected, 'Unlisted or missing files in transferred package: ' + str(sorted(actual ^ expected)[:8]))
    need(sum(int(r['bytes']) for r in records) == info['payload_bytes'], 'Payload size mismatch')
    print('VERIFIED', root, 'FILES', len(records), flush=True)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, default=DEFAULT_RUN)
    parser.add_argument('--project-root', type=Path, default=DEFAULT_PROJECT)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--dry-run', action='store_true', help='Validate selected files and show size; do not copy')
    parser.add_argument('--no-source', action='store_true', help='Collect results without a matching execution-source snapshot')
    parser.add_argument('--verify', type=Path, help='Verify a transferred collection instead of collecting')
    args = parser.parse_args(argv)
    previous = csv.field_size_limit()
    csv.field_size_limit(max(previous, 2**31-1))
    try:
        if args.verify:
            return verify_bundle(args.verify.expanduser().resolve())
        root = args.run.expanduser().resolve()
        project = None if args.no_source else args.project_root.expanduser().resolve()
        out = (args.output.expanduser().resolve() if args.output else
               default_review_output(root))
        need(root.is_dir(), 'Run directory not found: ' + str(root))
        need(not out.is_relative_to(root) and not root.is_relative_to(out),
             'Collection must be outside the experiment directory')
        if project:
            need(not out.is_relative_to(project) and not project.is_relative_to(out),
                 'Collection must be outside the code directory')
        need(not out.exists(), 'Collection output already exists: ' + str(out))
        print('READ ONLY RUN', root, flush=True)
        print('NO DOWNLOADS / NO ANALYSES / NO ZIP', flush=True)
        with read_locks(root):
            plan = build_plan(root, project)
            size = sum(source.stat().st_size for source, _ in plan.items.values())
            print(f'PLAN {len(plan.items)} files; {size / 1024**2:.1f} MiB', flush=True)
            print('COVERAGE 20 sections / 105 table blocks / 8 main figures / 12 supplementary figures', flush=True)
            for message in plan.warnings:
                print('NOTE', message, flush=True)
            if args.dry_run:
                return 0
            copy_bundle(plan, out)
        verify_bundle(out)
        print('REVIEW_DIRECTORY', out, flush=True)
        print('Transfer this directory only. Manuscript agreement remains to be checked.', flush=True)
        return 0
    finally:
        csv.field_size_limit(previous)


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        print('STOP:', exc, file=sys.stderr)
        raise SystemExit(2)
