#!/usr/bin/env python3
"""论文1公开数据分析入口；Linux运行，可核验复用本流程已经准备好的计数。"""
import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import time
import traceback
import zipfile
from pipeline.io import command, dump, require, sha

ROOT = Path(__file__).resolve().parent
VERSION = '1.2.2'


def fingerprint(config):
    """记录代码与参数，防止恢复时将两种分析混入同一个工作目录。"""
    digest = hashlib.sha256(json.dumps(config, sort_keys=True).encode())
    paths = [ROOT/'run.py', ROOT/'run_dream_repair.py']
    for folder in ('pipeline','analysis','tests'):
        paths.extend(p for p in (ROOT/folder).glob('*') if p.suffix in ('.py','.R'))
    for p in sorted(paths):
        digest.update(p.relative_to(ROOT).as_posix().encode())
        digest.update(p.read_bytes())
    return digest.hexdigest()


def dependencies():
    """先验证Python模块可导入，避免下载数百GB后才发现本机依赖缺失。"""
    for name in ('numpy', 'pandas', 'scipy', 'h5py', 'openpyxl', 'matplotlib'):
        module = importlib.import_module(name)
        print(name, module.__version__, flush=True)


def rrun(script, args, work, rscript):
    """在独立R进程运行一个阶段，并保存其完整报错。"""
    command([rscript, '--vanilla', ROOT/'analysis'/script, *args], work/'logs'/(script+'.log'))


def stage(name, function, outputs, work, resume):
    """只复用当前工作目录内、摘要一致且完成记录齐全的阶段。"""
    receipt = work/'stages'/(name+'.json')
    if receipt.exists():
        require(resume, f'阶段{name}已有完成记录；续跑请加--resume或使用新work目录')
        previous = json.loads(receipt.read_text())
        require(all((work/p).is_file() and sha(work/p)==digest for p,digest in previous['outputs'].items()),
                f'阶段{name}输出有变动；请在新目录重算')
        print('已核验，跳过', name, flush=True)
        return
    print('开始', name, flush=True)
    start = time.monotonic()
    function()
    paths = []
    for value in outputs:
        p = work/value
        require(p.exists(), f'阶段{name}没有生成{p}')
        paths.extend(sorted(x for x in p.rglob('*') if x.is_file()) if p.is_dir() else [p])
    require(paths, '阶段输出为空')
    dump(receipt, dict(stage=name, seconds=time.monotonic()-start,
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
    blocks, rows, _ = seaad.prepare_seaad(cache, folder)
    seaad.combine(blocks, rows, folder, 'human')
    quality_join(work/'raw', work/'prepared')


def prepare_mouse(work, config, library, rscript):
    """从六份GEO样本表和发表S2构建身份，再由R读取空间counts聚合。"""
    from pipeline.sources import mouse_manifest, unwrap_mouse
    mouse_manifest(work/'raw', work/'prepared')
    obj = unwrap_mouse(work/'raw')
    (work/'prepared'/'mouse_object_path.txt').write_text(str(obj)+'\n', encoding='utf-8')
    rrun('prepare_mouse.R', [work, library, ROOT], work, rscript)


def simulate(work, config, library, rscript):
    """从本次动物清单生成同场比较设计，再运行可断点续跑的R模拟。"""
    from pipeline.simulation_design import prepare
    prepare(work, config)
    rrun('simulation.R', [work, library, ROOT], work, rscript)


def bundle(work):
    """打包小型图表与状态，排除下载、聚合计数和运行缓存。"""
    path = work/'review_bundle.zip'
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as archive:
        for folder in ('figures', 'logs', 'stages', 'repair_details'):
            for p in (work/folder).rglob('*'):
                if p.is_file():
                    archive.write(p,p.relative_to(work))
        for name in ('status.json','run_config.json','software.json','input_reuse.json',
                     'repair_plan.json','repair_input_receipt.json','repair_changes.json','data_plan.json'):
            if (work/name).exists():
                archive.write(work/name,name)
        for p in (work/'results').rglob('*'):
            if p.is_file() and p.name != 'simulation_replicates.csv.gz':
                archive.write(p,p.relative_to(work))
        for name in ('human/metadata.csv','human/donor_inclusion.csv','human/source_object_inventory.csv',
                     'mouse_source_eligibility.csv','mouse_cage_labels.csv','mouse_inclusion.csv','mouse_units.csv','discovery_metadata.csv'):
            p=work/'prepared'/name
            if p.is_file():
                archive.write(p,p.relative_to(work))
        for p in (work/'prepared'/'simulation').glob('*.csv'):
            archive.write(p,p.relative_to(work))
        for p in (work/'raw').rglob('*.download.json'):
            archive.write(p,p.relative_to(work))
        # 附上实际执行的文本代码，下一次核对不再只能推测服务器用了哪个版本。
        code = [ROOT/'run.py', ROOT/'run_dream_repair.py', ROOT/'config.json']
        for folder in ('pipeline', 'analysis', 'tests'):
            code.extend(p for p in (ROOT/folder).glob('*') if p.suffix in ('.py', '.R'))
        for p in code:
            archive.write(p, 'executed_code/' + p.relative_to(ROOT).as_posix())
        archive.writestr('executed_code/sha256.json', json.dumps(
            {p.relative_to(ROOT).as_posix(): sha(p) for p in code}, indent=2))
    print('审阅包：', path, flush=True)


def main():
    """命令行调度。安装必须显式请求；任何阶段失败返回非零并生成错误日志。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', choices=['all','reanalyse','check','install-r','download','prepare','models','simulation','calibration','figures','audit'], default='all')
    parser.add_argument('--work', type=Path, default=ROOT/'work_public')
    parser.add_argument('--config', type=Path, default=ROOT/'config.json')
    parser.add_argument('--rscript', default=shutil.which('Rscript') or 'Rscript')
    parser.add_argument('--library', type=Path, default=ROOT/'R_library')
    parser.add_argument('--workers', type=int)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--reuse-work', '--data-work', dest='reuse_work', type=Path,
                        help='已有数据工作目录：all优先复用合格计数和原始缓存，仅下载缺失文件；reanalyse要求完整聚合输入')
    args = parser.parse_args()
    print('运行脚本版本：', VERSION, flush=True)
    require(platform.system()=='Linux', '本程序只在Linux运行；本地单元测试可用python -m unittest')
    config = json.loads(args.config.read_text(encoding='utf-8'))
    if args.workers is not None:
        config['workers'] = args.workers
    require(config['workers']>=1 and config['calibration_repeats']>=100 and config['simulation_repeats']>=100, '并行数或重复数不合法')
    require(config.get('analysis_version') == VERSION, '请使用本版config.json')
    require(isinstance(config.get('ebayes_robust'), bool) and isinstance(config.get('dream_ebayes_robust'), bool),
            '必须明确填写经验贝叶斯设置')
    require(args.task != 'reanalyse' or args.reuse_work is not None, 'reanalyse须指定--data-work')
    require(args.reuse_work is None or args.task in ('all','reanalyse'), '--data-work仅用于all或reanalyse')
    work = args.work.expanduser().resolve()
    if args.reuse_work:
        source = args.reuse_work.expanduser().resolve()
        require(source != work and not work.is_relative_to(source) and not source.is_relative_to(work),
                '新结果目录必须与已有数据目录分开，不能相同或相互嵌套')
    library = args.library.expanduser().resolve()
    signature = fingerprint(config)
    if work.exists() and any(work.iterdir()):
        require((work/'owner.json').exists(), 'work目录非空且不是本程序创建；请指定全新目录')
    work.mkdir(parents=True,exist_ok=True)
    for folder in ('raw','prepared','results','logs','stages'):
        (work/folder).mkdir(exist_ok=True)
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
        if (work/'owner.json').exists():
            previous = json.loads((work/'owner.json').read_text())['signature']
            require(previous == signature, '代码或参数已变；请用--task reanalyse --reuse-work 原目录 --work 新目录，不要删除旧记录')
    except Exception:
        lock.close()
        raise
    try:
        dump(work/'owner.json',dict(signature=signature,source='official_public_downloads'))
        dump(work/'run_config.json',config)
        os.environ['PYTHONPATH'] = str(ROOT) + os.pathsep + os.environ.get('PYTHONPATH','')
        dump(work/'status.json',dict(status='running',task=args.task,pid=os.getpid()))
        dependencies()
        dump(work/'software.json',dict(python=sys.version,platform=platform.platform()))
        if args.task=='install-r':
            rrun('preflight.R',['install',library,ROOT],work,args.rscript)
            rrun('preflight.R',['check',library,ROOT,work/'logs'],work,args.rscript)
            dump(work/'status.json',dict(status='dependencies_passed',task=args.task,analysis_run=False))
            return 0
        command([sys.executable,'-m','unittest','discover','-s',ROOT/'tests','-v'],work/'logs'/'tests.log')
        rrun('preflight.R',['check',library,ROOT,work/'logs'],work,args.rscript)
        if args.task=='check':
            dump(work/'status.json',dict(status='checks_passed',task=args.task,analysis_run=False))
            return 0
        reused = set()
        if args.task == 'reanalyse':
            from pipeline.reuse import import_inputs
            import_inputs(args.reuse_work, work, config)
        elif args.task == 'all' and args.reuse_work is not None:
            from pipeline.cache import prepare_existing_data
            reused = prepare_existing_data(args.reuse_work, work, config)
        from pipeline.sources import fetch_small, prepare_panel, map_orthologs
        plan = [
            ('download', lambda: fetch_small(work/'raw'), ['raw']),
            ('prepare_panel', lambda: prepare_panel(work/'raw',work/'prepared'), ['prepared/discovery_expression.csv.gz','prepared/discovery_metadata.csv','prepared/probe_mapping.csv']),
            ('panel', lambda: rrun('panel.R',[work,library,ROOT],work,args.rscript), ['results/panel.csv','results/panel_discovery.csv.gz']),
            ('orthologs', lambda: map_orthologs(work/'raw',work/'results'), ['results/ortholog_panel.csv']),
            ('prepare_human', lambda: prepare_human(work,config), ['prepared/human']),
            ('prepare_mouse', lambda: prepare_mouse(work,config,library,args.rscript), ['prepared/mouse_units.csv','prepared/mouse_whole_counts.csv.gz','prepared/mouse_cortex_counts.csv.gz']),
            ('human', lambda: rrun('models.R',['human',work,library,ROOT],work,args.rscript), ['results/human_results.csv.gz','results/human_paired_changes.csv.gz','results/correction_scope.csv','results/human_designs.csv','results/human_gene_membership.csv','results/human_correction_scope_summary.csv']),
            ('mouse', lambda: rrun('models.R',['mouse',work,library,ROOT],work,args.rscript),
             ['results/mouse_results.csv.gz','results/mouse_precision.csv.gz','results/mouse_variance.csv.gz',
              'results/mouse_designs.csv','results/mouse_background.csv','results/mouse_SE_implementation_audit.csv.gz',
              'results/mouse_background_fits.csv.gz','results/mouse_fit_failures.csv',
              'results/mouse_model_status.csv','results/mouse_model_messages.csv','results/mouse_dream_diagnostics.csv',
              'results/mouse_dream_retry_attempts.csv','results/mouse_dream_background_impact.csv.gz']),
        ]
        from pipeline.calibration import run as calibrate
        from pipeline.figures import generate
        from pipeline.audit import assess
        plan += [
            ('simulation', lambda: simulate(work,config,library,args.rscript),
             ['results/simulation_summary.csv','results/simulation_replicates.csv.gz','results/simulation_layouts.csv',
              'results/simulation_scenarios.csv','results/simulation_failures.csv','prepared/simulation']),
            ('calibration', lambda: calibrate(work/'prepared'/'human'/'metadata.csv',work/'results'/'calibration',config['calibration_repeats'],config['seed']), ['results/calibration']),
            ('figures', lambda: generate(work), ['figures']),
            ('audit', lambda: assess(work), ['results/result_audit.json','results/manuscript_comparison.csv'])]
        selected = {'download':['download'], 'prepare':['prepare_panel','panel','orthologs','prepare_human','prepare_mouse'],
                    'models':['human','mouse'], 'simulation':['simulation'], 'calibration':['calibration'], 'figures':['figures'],
                    'audit':['audit'], 'reanalyse':['human','mouse','simulation','calibration','figures','audit']}
        for name, function, outputs in plan:
            if args.task=='all' or name in selected.get(args.task,[]):
                if name in reused or (name=='prepare_panel' and 'panel' in reused) or (
                    name=='download' and reused=={'panel','orthologs','prepare_human','prepare_mouse'}):
                    print('已有合格输入，跳过下载或数据准备：',name,flush=True)
                    continue
                stage(name,function,outputs,work,args.resume)
        audit = json.loads((work/'results'/'result_audit.json').read_text()) if args.task in ('all','reanalyse','audit') else {}
        dump(work/'status.json',dict(status=audit.get('status','requested_stages_completed'),task=args.task,version=VERSION,
             audit=audit,
             historical_manuscript_reproduction_certified=False,
             note='本次公开数据重分析完成。与历史论文逐项核对是后续独立步骤。'))
        bundle(work)
        if audit.get('status') in ('completed_with_fit_failures','completed_with_unresolved_warnings'):
            print('分析已结束，仍需审阅拟合或收敛明细：', audit, flush=True)
            print('退出码3表示结果不宜直接定稿，不是下载或进程崩溃。', flush=True)
            return 3
        return 0
    except Exception as exc:
        (work/'logs'/'error.txt').write_text(traceback.format_exc(),encoding='utf-8')
        dump(work/'status.json',dict(status='failed',task=args.task,reason=str(exc)))
        bundle(work)
        print('停止：',exc,flush=True)
        return 2
    finally:
        lock.close()


if __name__=='__main__':
    sys.exit(main())
