#!/usr/bin/env python3
"""Public expression inputs and primary fitted-model analyses."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import time
import traceback
from pipeline.io import command, dump, require, sha
from pipeline.configuration import VERSION, load_config

ROOT = Path(__file__).resolve().parents[1]


def fingerprint(config):
    """记录代码与参数，防止恢复时将两种分析混入同一个工作目录。"""
    digest = hashlib.sha256(json.dumps(config, sort_keys=True).encode())
    paths = [ROOT/'run.py']
    for folder in ('pipeline','analysis'):
        paths.extend(p for p in (ROOT/folder).rglob('*') if p.suffix in ('.py','.R'))
    for p in sorted(paths):
        digest.update(p.relative_to(ROOT).as_posix().encode())
        digest.update(p.read_bytes())
    return digest.hexdigest()


def dependencies():
    """先验证Python模块可导入，避免下载数百GB后才发现本机依赖缺失。"""
    from pipeline.environment import python_snapshot
    record = python_snapshot()
    for name, version in record['packages'].items():
        print(name, version, flush=True)
    return record


def rrun(script, args, work, rscript):
    """在独立R进程运行一个阶段，并保存其完整报错。"""
    command([rscript, '--vanilla', ROOT/'analysis'/script, *args], work/'logs'/(script+'.log'))


def stage(name, function, outputs, work, resume):
    """只复用当前工作目录内、摘要一致且完成记录齐全的阶段。"""
    receipt = work/'stages'/(name+'.json')
    from pipeline.lifecycle import assert_code
    identity = json.loads((work/'fresh_run.json').read_text(encoding='utf-8'))['code']
    assert_code(ROOT, identity)
    if receipt.exists():
        require(resume, f'阶段{name}已有完成记录；续跑请加--resume或使用新work目录')
        previous = json.loads(receipt.read_text())
        require(previous.get('signature') == json.loads((work/'owner.json').read_text())['signature'],
                'Stage receipt is not from this raw-input experiment')
        records = previous.get('outputs')
        require(isinstance(records, dict) and records, f'阶段{name}完成记录没有输出文件')
        require(all(not Path(p).is_absolute() and '..' not in Path(p).parts
                    and (work/p).resolve().is_relative_to(work.resolve()) for p in records),
                'Stage receipt points outside the current experiment')
        require(all((work/p).is_file() and sha(work/p)==digest for p,digest in records.items()),
                f'阶段{name}输出有变动；请在新目录重算')
        print('已核验，跳过', name, flush=True)
        return
    print('开始', name, flush=True)
    start = time.monotonic()
    function()
    assert_code(ROOT, identity)
    paths = []
    for value in outputs:
        p = work/value
        require(p.exists(), f'阶段{name}没有生成{p}')
        paths.extend(sorted(x for x in p.rglob('*') if x.is_file()) if p.is_dir() else [p])
    require(paths, '阶段输出为空')
    dump(receipt, dict(stage=name, signature=json.loads((work/'owner.json').read_text())['signature'],
                       seconds=time.monotonic()-start,
                       outputs={p.relative_to(work).as_posix():sha(p) for p in paths}))


def prepare_human(work, config):
    """下载指定SEA-AD发行版全部DFC供体，按亚类聚合原始UMI并连接质量信息。"""
    from pipeline import seaad
    from pipeline.io import download
    from pipeline.sources import quality_join
    seaad.RELEASE = config['seaad_release']
    seaad.EXPECTED_DONORS = config['seaad_expected_donors']
    seaad.download = download
    folder = work/'prepared'/'human'
    folder.mkdir(exist_ok=True)
    cache = work/'raw'/'SEAAD'
    cache.mkdir(exist_ok=True)
    run_id = json.loads((work/'fresh_run.json').read_text())['run_id']
    blocks, rows, _ = seaad.prepare_seaad(cache, folder, run_id=run_id)
    seaad.combine(blocks, rows, folder, 'human')
    quality_join(work/'raw', work/'prepared')


def prepare_mouse(work, config, library, rscript):
    """从六份GEO样本表和发表S2构建身份，再由R读取空间counts聚合。"""
    from pipeline.sources import mouse_manifest, unwrap_mouse
    mouse_manifest(work/'raw', work/'prepared')
    obj = unwrap_mouse(work/'raw')
    (work/'prepared'/'mouse_object_path.txt').write_text(str(obj)+'\n', encoding='utf-8')
    rrun('prepare_mouse.R', [work, library, ROOT], work, rscript)


def main():
    """命令行调度。安装必须显式请求；任何阶段失败返回非零并生成错误日志。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', choices=['all','check','install-r','download','prepare','models','calibration'], default='all')
    parser.add_argument('--work', type=Path, default=ROOT/'work_public')
    parser.add_argument('--config', type=Path, default=ROOT/'config.json')
    parser.add_argument('--rscript', default=shutil.which('Rscript') or 'Rscript')
    parser.add_argument('--library', type=Path, default=ROOT/'R_library')
    parser.add_argument('--workers', type=int)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--raw-cache', type=Path, action='append', default=[],
                        help='Reuse only registered public downloads, never prepared inputs or results')
    args = parser.parse_args()
    print('运行脚本版本：', VERSION, flush=True)
    require(platform.system()=='Linux', 'Scientific execution requires Linux; Python unit tests are cross-platform')
    config = load_config(args.config, args.workers)
    work = args.work.expanduser().resolve()
    library = args.library.expanduser().resolve()
    work.mkdir(parents=True,exist_ok=True)
    library.mkdir(parents=True,exist_ok=True)
    # Linux锁在进程退出后由内核释放，不靠陈旧PID文件判断是否仍在运行。
    import fcntl
    lock = (work/'run.lock').open('w')
    try:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError('同一work目录已有运行进程')
    # 代码身份检查失败时，不覆盖已有状态；新版本使用不同输出目录。
    try:
        from pipeline.lifecycle import claim, parent_run, code_identity
        run_id = claim(work, 'core', args.resume, parent_run(), code_identity(ROOT))
        signature = fingerprint(dict(run_id=run_id, settings={k:v for k,v in config.items() if k != 'workers'}))
        if (work/'owner.json').exists():
            previous = json.loads((work/'owner.json').read_text())['signature']
            require(previous == signature, 'Scientific code or settings changed; use a NEW output and --raw-cache only')
    except Exception:
        lock.close()
        raise
    try:
        for folder in ('raw','prepared','results','logs','stages'):
            (work/folder).mkdir(exist_ok=True)
        dump(work/'owner.json',dict(signature=signature,source='official_public_downloads'))
        dump(work/'run_config.json',config)
        os.environ['PYTHONPATH'] = str(ROOT) + os.pathsep + os.environ.get('PYTHONPATH','')
        dump(work/'status.json',dict(status='running',task=args.task,pid=os.getpid()))
        python_environment = dependencies()
        if args.task=='install-r':
            rrun('preflight.R',['install',library,ROOT],work,args.rscript)
            rrun('preflight.R',['check',library,ROOT,work/'logs'],work,args.rscript)
            dump(work/'status.json',dict(status='dependencies_passed',task=args.task,analysis_run=False))
            return 0
        from pipeline.environment import r_snapshot, freeze_environment
        freeze_environment(work/'software.json', {'python': python_environment,
                           'R': r_snapshot(args.rscript, [library])})
        command([sys.executable,'-m','unittest','discover','-s',ROOT/'tests','-v'],work/'logs'/'tests.log')
        rrun('preflight.R',['check',library,ROOT,work/'logs'],work,args.rscript)
        if args.task=='check':
            dump(work/'status.json',dict(status='checks_passed',task=args.task,analysis_run=False))
            return 0
        from pipeline.cache import import_raw_sources
        import_raw_sources(args.raw_cache, work, config)
        from pipeline.sources import fetch_small, prepare_panel, map_orthologs
        plan = [
            ('download', lambda: fetch_small(work/'raw'), ['raw']),
            ('prepare_panel', lambda: prepare_panel(work/'raw',work/'prepared'), ['prepared/discovery_expression.csv.gz','prepared/discovery_metadata.csv','prepared/probe_mapping.csv']),
            ('panel', lambda: rrun('panel.R',[work,library,ROOT],work,args.rscript), ['results/panel.csv','results/panel_discovery.csv.gz']),
            ('orthologs', lambda: map_orthologs(work/'raw',work/'results'), ['results/ortholog_panel.csv']),
            ('prepare_human', lambda: prepare_human(work,config), ['prepared/human']),
            ('prepare_mouse', lambda: prepare_mouse(work,config,library,args.rscript), ['prepared/mouse_units.csv','prepared/mouse_whole_counts.csv.gz','prepared/mouse_cortex_counts.csv.gz']),
            ('human', lambda: rrun('models.R',['human',work,library,ROOT],work,args.rscript), ['results/human_results.csv.gz','results/human_paired_changes.csv.gz','results/human_precision.csv.gz','results/correction_scope.csv','results/human_designs.csv','results/human_gene_membership.csv','results/human_correction_scope_summary.csv']),
            ('mouse', lambda: rrun('models.R',['mouse',work,library,ROOT],work,args.rscript),
             ['results/mouse_results.csv.gz','results/mouse_precision.csv.gz','results/mouse_variance.csv.gz',
              'results/mouse_designs.csv','results/mouse_background.csv','results/mouse_SE_decomposition_check.csv.gz',
              'results/mouse_background_fits.csv.gz','results/mouse_fit_failures.csv',
              'results/mouse_model_status.csv','results/mouse_model_messages.csv','results/mouse_dream_diagnostics.csv',
              'results/mouse_dream_retry_attempts.csv','results/mouse_dream_background_impact.csv.gz']),
        ]
        from pipeline.calibration import run as calibrate
        plan.append(('calibration', lambda: calibrate(work/'prepared'/'human'/'metadata.csv',
                      work/'results'/'calibration',config['calibration_repeats'],config['seed']),
                      ['results/calibration']))
        selected = {'download':['download'],
                    'prepare':['download','prepare_panel','panel','orthologs','prepare_human','prepare_mouse'],
                    'models':['human','mouse'], 'calibration':['calibration']}
        for name, function, outputs in plan:
            if args.task=='all' or name in selected.get(args.task,[]):
                stage(name,function,outputs,work,args.resume)
        dump(work/'status.json',dict(status='completed',task=args.task,version=VERSION))
        return 0
    except Exception as exc:
        (work/'logs'/'error.txt').write_text(traceback.format_exc(),encoding='utf-8')
        dump(work/'status.json',dict(status='failed',task=args.task,reason=str(exc)))
        print('停止：',exc,flush=True)
        return 2
    finally:
        lock.close()


if __name__=='__main__':
    sys.exit(main())
