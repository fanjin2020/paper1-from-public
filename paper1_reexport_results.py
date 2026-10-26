#!/usr/bin/env python3
"""Re-export one completed run into a new review directory, without model refits."""
import argparse
import csv
import hashlib
from datetime import datetime, timezone
import importlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

PROJECT = Path(__file__).resolve().parent
DEFAULT_RUN = PROJECT.parent / 'paper1_complete_from_public_20261022_v2'
PACKAGES = ('numpy', 'pandas', 'scipy', 'matplotlib')


def environment(record):
    errors = []
    if platform.python_version() != record.get('python'):
        errors.append('Python: expected '+str(record.get('python'))+', found '+platform.python_version())
    versions = {}
    for name in PACKAGES:
        try:
            value = str(importlib.import_module(name).__version__)
            versions[name] = value
            if value != record.get('packages', {}).get(name):
                errors.append(name+': expected '+str(record.get('packages', {}).get(name))+', found '+value)
        except Exception as exc:
            errors.append(name+': '+str(exc))
    if errors:
        raise ValueError('Reporting environment differs from the completed run; nothing installed.\n'+'\n'.join(errors))
    return dict(python=platform.python_version(), executable=sys.executable, packages=versions,
                unused='R, openpyxl, h5py, scrublet and raw-data preparation packages')


def safe_output(root, out, project=PROJECT):
    for other in (root, project):
        if out.is_relative_to(other) or other.is_relative_to(out):
            raise ValueError('Export output must be outside the original run and code directories')
    if out.exists():
        raise ValueError('Output already exists. Use a new directory; originals are never overwritten: '+str(out))


def input_plan(root, collector):
    if not (root/'COLLECTION_COMPLETE.json').is_file():
        return collector.build_plan(root, project=None), False
    # Portable collections retain server paths in immutable receipts. Accept
    # relocation only after verifying the entire checksummed file inventory.
    collector.verify_bundle(root)
    plan=collector.Plan(root)
    collector.validate_completion(plan)
    for row in collector.csv_records(root/'file_manifest.csv'):
        name=row['file']
        if name.startswith(('core/','matched/','benchmark/','supplementary/','maintenance/','source_snapshot/')) or name in (
            'state.json','software.json','fresh_run.json','paper_coverage.json'):
            plan.add(name,row['category'])
            plan.expect(collector.inside(root,name),row['sha256'],int(row['bytes']))
    return plan, True


def export(root, out, runtime):
    from pipeline import result_collection as C
    from pipeline import reporting, figures, supplementary_figures, publication_tables
    from pipeline.lifecycle import code_identity
    safe_output(root, out)
    def identity():
        return hashlib.sha256((code_identity(PROJECT)+C.sha256(Path(__file__).resolve())).encode()).hexdigest()
    export_code = identity()
    with C.read_locks(root):
        print('CHECK COMPLETED INPUTS AND RECEIPTS', flush=True)
        plan, portable = input_plan(root, C)
        # Completion receipts describe the scientific run. Updated export code
        # has its own identity and must not rewrite those historical identities.
        plan.items = {name:item for name,item in plan.items.items()
                      if not name.startswith(('reporting/', 'figures/'))}
        out.mkdir(parents=True)
        C.write_json(out/'COLLECTION_IN_PROGRESS.json', dict(status='EXPORTING', source_run=str(root)))
        try:
            manifest = []
            for i,(name,(source,category)) in enumerate(sorted(plan.items.items()),1):
                target=C.inside(out,name);target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(source,target)
                expected,size=plan.expected[str(source)]
                C.need(target.stat().st_size==size and C.sha256(target)==expected, 'Copy checksum mismatch: '+name)
                if i%100==0:print('COPY COMPLETED RESULTS',i,'/',len(plan.items),flush=True)
            supplements={name:root/'supplementary'/name/Path(C.read_json(root/'supplementary'/(name+'.complete.json'))['directory']).name
                         for name in C.GROUPS}
            for name,path in supplements.items():
                recorded=Path(C.read_json(root/'supplementary'/(name+'.complete.json'))['directory']).resolve()
                C.need((portable or recorded==path.resolve()) and path.is_dir(), 'Supplementary attempt mismatch: '+name)
            report=out/'reporting';report.mkdir()
            print('REBUILD 20 WORKSHEETS / 112 EXPLICIT BLOCKS / XLSX',flush=True)
            index=reporting.generate(root,report,supplements)
            plan.checks['source_table_blocks']=C.validate_block_mapping(list(C.csv_records(root/'reporting/block_sources.csv')))
            plan.checks['table_blocks']=len(publication_tables.EXPECTED_BLOCKS)
            print('DRAW 12 SUPPLEMENTARY FIGURES',flush=True)
            plots=supplementary_figures.generate(report)
            print('DRAW 8 MAIN FIGURES',flush=True)
            figures.generate(root/'core',root/'matched',root/'benchmark',out/'figures')
            C.write_json(report/'status.json',dict(status='RESULT_EXPORTS_COMPLETE',run_id=plan.checks['run_id'],
                files=len(index),sections=20,workbook_blocks_mapped=len(publication_tables.EXPECTED_BLOCKS),supplementary_figures=len(plots),
                workbook='Supplementary_Tables.xlsx',workbook_sha256=C.sha256(report/'Supplementary_Tables.xlsx'),
                reporting_version=publication_tables.VERSION,model_refits=False,manuscript_agreement='NOT_ASSESSED'))
            C.need(identity()==export_code, 'Export source changed while running')
            for folder in ('pipeline','analysis','tests','inputs','docs'):
                for source in sorted((PROJECT/folder).rglob('*')):
                    if source.suffix in ('.py','.R','.json','.csv','.md'):
                        target=out/'export_source_snapshot'/source.relative_to(PROJECT)
                        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
            for name in ('paper1_reexport_results.py','run.py','config.json','requirements.txt','requirements-supplementary.txt',
                         'reproduce_from_public.py','README.md','README_CN.md'):
                source=PROJECT/name
                if source.is_file():shutil.copyfile(source,out/'export_source_snapshot'/name)
            # Preserve unchanged scientific results, including genuine differences
            # from the manuscript; never substitute numbers from an earlier run.
            C.write_json(out/'export_provenance.json',dict(source_run=str(root),run_id=plan.checks['run_id'],
                scientific_source_identity=plan.checks['source_identity'],export_source_identity=export_code,
                reporting_version=publication_tables.VERSION,environment=runtime,model_refits=False,downloads=False,
                original_run_modified=False,input_is_verified_portable_collection=portable,
                statistical_differences_not_repaired_by_export=True,
                source_snapshot_role='EXPORT_CODE; scientific execution identity remains in original completion records'))
            C.write_json(out/'collection_details.json',dict(source_run=str(root),checks=plan.checks,
                warnings=['This export corrects reporting definitions, not scientific fit differences.',
                          'Execution code identity and export code identity are intentionally separate.'],
                excluded=plan.excluded))
            (out/'README_REVIEW.md').write_text(
                '# Paper1 re-exported verification materials\n\n'
                '20 worksheets, 112 explicit block sources, 8 main and 12 supplementary figures.\n'
                'reporting/Supplementary_Tables.xlsx contains all 112 sections with an S01 catalog.\n'
                'No fitting, simulations or downloads were performed. The workbook was assembled from computed results.\n'
                'Original run files remain unchanged. Scientific tables in core/matched/benchmark/supplementary are byte-for-byte copies.\n'
                'reporting contains derived publication tables. Consult block_sources.csv for each block.\n'
                'S4.2 summarizes all finite clinical log pairs; selection is counted separately.\n'
                'S9 interaction signs are not disease-direction comparisons; original source copies retain their recorded labels.\n'
                'S19.4 distinguishes actual five-model and three-model layout membership.\n'
                'S15.1 reports absolute coefficient-join differences and a separate signed column.\n'
                'Native two-factor precision and shared-weight three-factor precision are named separately.\n'
                'ICC uses dream/base SE ratios with .01/.10/.20 bin boundaries.\n'
                'Real GSE222007/GSE167490/mixture/S18 differences from the manuscript remain for review.\n'
                'export_source_snapshot is the new exporter, not the original fitting-code snapshot.\n'
                'No claim of complete manuscript agreement is made. Transfer this entire directory to a short parent path.\n',encoding='utf-8')
            shutil.copyfile(PROJECT/'pipeline/result_collection.py',out/'paper1_collect_verification.py')
            for path in sorted(out.rglob('*')):
                if not path.is_file() or path.name=='COLLECTION_IN_PROGRESS.json':continue
                rows,columns=C.table_info(path)
                manifest.append(dict(file=path.relative_to(out).as_posix(),category='reexport',bytes=path.stat().st_size,
                    sha256=C.sha256(path),rows=rows,columns=columns))
            C.write_csv(out/'file_manifest.csv',manifest,['file','category','bytes','sha256','rows','columns'])
            # Detect source changes during copy/summary without reading raw counts.
            for name,(source,_) in plan.items.items():
                digest,size=plan.expected[str(source)]
                C.need(source.stat().st_size==size and C.sha256(source)==digest, 'Source changed during export: '+name)
            (out/'COLLECTION_IN_PROGRESS.json').unlink()
            C.write_json(out/'COLLECTION_COMPLETE.json',dict(status='VERIFICATION_PACKAGE_COMPLETE',
                run_id=plan.checks['run_id'],created_utc=datetime.now(timezone.utc).isoformat(),files=len(manifest),
                payload_bytes=sum(r['bytes'] for r in manifest),manifest_sha256=C.sha256(out/'file_manifest.csv'),
                manuscript_agreement='NOT_ASSESSED',model_refits=False))
            C.verify_bundle(out)
        except BaseException as exc:
            C.write_json(out/'COLLECTION_FAILED.json',dict(status='INCOMPLETE',reason=str(exc)))
            raise
    print('REVIEW_DIRECTORY',out,flush=True)
    print('Transfer this directory only. No ZIP and no raw counts.',flush=True)
    return 0


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,default=DEFAULT_RUN)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--python',type=Path,help='Existing interpreter matching the recorded reporting environment')
    parser.add_argument('--verify',type=Path,help='Verify the exported directory without scientific packages')
    parser.add_argument('--execute',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args(argv)
    csv.field_size_limit(max(csv.field_size_limit(),2**31-1))
    if args.verify:
        from pipeline.result_collection import verify_bundle
        return verify_bundle(args.verify.resolve())
    root=args.run.expanduser().resolve()
    from pipeline.result_collection import default_review_output
    out=args.output.expanduser().resolve() if args.output else default_review_output(root)
    safe_output(root,out)
    recorded=json.loads((root/'software.json').read_text(encoding='utf-8'))['python']
    if not args.execute:
        interpreter=str(args.python.expanduser() if args.python else recorded['executable'])
        if not Path(interpreter).is_file():raise ValueError('Recorded Python unavailable; supply --python for the existing matching environment: '+interpreter)
        env=os.environ.copy();env['PYTHONPATH']=str(PROJECT)
        for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
        print('REUSE RECORDED PYTHON',interpreter,'\nNO DOWNLOADS / NO MODEL REFITS / NO ZIP',flush=True)
        return subprocess.call([interpreter,'-u',str(Path(__file__).resolve()),'--execute','--run',str(root),'--output',str(out)],cwd=PROJECT,env=env)
    runtime=environment(recorded)
    print('REPORTING ENVIRONMENT',runtime,flush=True)
    return export(root,out,runtime)


if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,OSError,KeyError) as error:
        print('STOP:',error,file=sys.stderr);raise SystemExit(2)
