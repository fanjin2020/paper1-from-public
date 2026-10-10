#!/usr/bin/env python3
"""核验v1.2.0输出，只复核出现数值警告的dream模型，不下载或重跑其它分析。"""
import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import traceback

import numpy as np
import pandas as pd

from pipeline.io import command, dump, require, sha
from pipeline.reuse import checked_path
from pipeline.audit import assess
from run import ROOT, VERSION, bundle, fingerprint, rrun


def verify_release(root):
    """核对当前包的文件清单，避免把两个版本的修复文件混在一起。"""
    record = json.loads((root / 'release_manifest.json').read_text(encoding='utf-8'))
    require(record['version'] == VERSION, '发布清单版本不一致')
    for name, digest in record['files'].items():
        require(sha(checked_path(root, name)) == digest, f'代码文件与发布包不一致：{name}')


def verify_source(source):
    """按来源阶段收据核对全部结果及必需计数，不接受只有review_bundle的小包。"""
    cfg = json.loads((source / 'run_config.json').read_text(encoding='utf-8'))
    require(cfg.get('analysis_version') == '1.2.0', '定点修复入口只接受本次v1.2.0完整工作目录')
    status = json.loads((source / 'status.json').read_text(encoding='utf-8'))
    require(status.get('status') == 'completed_with_fit_failures', '来源并非本次带拟合失败的完成运行')
    require(status.get('audit', {}).get('integrity_checks') == 'passed', '来源尚未通过完整性核查')
    registry = {}
    stages = ('panel', 'orthologs', 'human', 'mouse', 'simulation', 'calibration', 'prepare_mouse')
    for stage in stages:
        receipt = source / 'stages' / f'{stage}.json'
        item = json.loads(receipt.read_text(encoding='utf-8'))
        require(item.get('stage') == stage, f'阶段名称不一致：{stage}')
        for name, digest in item['outputs'].items():
            require(name not in registry or registry[name] == digest, f'收据冲突：{name}')
            registry[name] = digest
    names = [p.relative_to(source).as_posix() for p in (source / 'results').rglob('*') if p.is_file()]
    names = [n for n in names if n not in ('results/result_audit.json', 'results/manuscript_comparison.csv')]
    names += ['prepared/mouse_units.csv', 'prepared/mouse_whole_counts.csv.gz',
              'prepared/mouse_cortex_counts.csv.gz']
    # sessionInfo不在旧mouse阶段的输出清单里，记录其摘要并在R中核对实际包版本。
    unregistered = {'results/human_R_session.txt', 'results/mouse_R_session.txt'}
    verified = {}
    for name in sorted(set(names)):
        p = checked_path(source, name)
        require(name in registry or name in unregistered, f'结果缺少来源收据：{name}')
        digest = sha(p)
        if name in registry:
            require(digest == registry[name], f'来源内容已改变：{name}')
        verified[name] = digest
    failures = pd.read_csv(source / 'results/mouse_fit_failures.csv')
    diagnostics = pd.read_csv(source / 'results/mouse_dream_diagnostics.csv')
    require(len(failures) > 0, '没有需要定点复核的失败记录')
    require(failures.model.isin(['dream_native', 'dream_shared']).all(), '存在非dream失败，不能用本入口替代完整排查')
    require(failures.status.eq('numerical_warning').all(), '存在不同类型的失败，需要单独核查')
    require(not diagnostics.gene_symbol.isna().any(), '存在未定位的警告，不能仅重试已知基因')
    require(set(zip(failures.spec, failures.model, failures.gene_symbol)) <=
            set(zip(diagnostics.spec, diagnostics.model, diagnostics.gene_symbol)), '失败记录与警告定位不一致')
    return cfg, verified, failures


def import_reviewed(source, work, verified):
    """复制结果；大计数使用只读访问的符号链接，绝不改写来源目录中的数据。"""
    for name, digest in verified.items():
        src, dst = source / name, work / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        if name.endswith('_counts.csv.gz'):
            dst.symlink_to(src)
        else:
            shutil.copy2(src, dst)
        require(sha(dst) == digest, f'复制后摘要不一致：{name}')
    dump(work / 'repair_input_receipt.json', dict(source=str(source), files=verified))


def bh_planned(values):
    """按计划位置计算BH，缺失值占分母但输出仍为缺失。"""
    p = np.asarray(values, dtype=float)
    ok = np.isfinite(p)
    require(((p[ok] >= 0) & (p[ok] <= 1)).all(), 'p值超出范围')
    filled = np.where(ok, p, 1.0)
    order = np.argsort(filled, kind='stable')
    ranked = np.minimum.accumulate((filled[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    q = np.empty(len(p)); q[order] = np.minimum(ranked, 1.0)
    q[~ok] = np.nan
    return q


def merge_repaired(work, targets):
    """替换完整的受影响模型，再对原计划面板共同重算BH；其它模型数值保持不变。"""
    r = work / 'results'
    bg = pd.read_csv(r / 'mouse_background_fits.csv.gz')
    replacement = pd.read_csv(r / 'dream_repair_replacement_rows.csv.gz')
    keys = ['spec', 'model', 'gene_symbol']
    require(not replacement.duplicated(keys).any(), '重拟合结果有重复位置')
    old = bg.set_index(keys)
    new = replacement.set_index(keys)
    require(new.index.isin(old.index).all(), '重拟合结果出现计划外基因')
    selected = pd.MultiIndex.from_frame(bg[['spec', 'model']]).isin(targets)
    require(len(new) == int(selected.sum()), '重拟合没有覆盖该模型的全部背景')
    require(set(new.index) == set(old.index[selected]), '背景位置不完整')
    changes = old.loc[new.index].add_suffix('_before').join(new.add_suffix('_after'))
    changes.to_csv(r / 'dream_repair_before_after.csv.gz', compression='gzip')
    old.loc[new.index, new.columns] = new
    bg = old.reset_index()
    bg.to_csv(r / 'mouse_background_fits.csv.gz', index=False, compression='gzip')
    bg.loc[~bg.status.eq('available')].to_csv(r / 'mouse_fit_failures.csv', index=False)
    summary = bg.groupby(['spec', 'model'], sort=False).agg(planned=('status','size'),
                    available=('status', lambda x: int(x.eq('available').sum()))).reset_index()
    summary['failed'] = summary.planned-summary.available
    summary.to_csv(r / 'mouse_model_status.csv', index=False)
    panel = pd.read_csv(r / 'mouse_results.csv.gz')
    numeric = ['beta','SE','df','p','low','high','status']
    lookup = new
    for i, row in panel.iterrows():
        key = (row.spec, row.model, row.mouse_gene)
        if key in lookup.index:
            panel.loc[i, numeric] = lookup.loc[key, numeric].to_numpy()
    panel['q_within_model'] = panel.groupby('model', sort=False).p.transform(bh_planned)
    for name, models in {
        'q_base_fixed_joint':['baseline','cage_fixed'],
        'q_three_model_joint':['baseline','cage_fixed','cage_correlation'],
        'q_dream_joint':['dream_shared','dream_native'],
    }.items():
        k = panel.model.isin(models)
        panel[name] = np.nan
        panel.loc[k,name] = bh_planned(panel.loc[k,'p'])
    panel.to_csv(r / 'mouse_results.csv.gz', index=False, compression='gzip')
    for original, incoming in [('mouse_dream_diagnostics.csv','dream_repair_diagnostics.csv'),
                               ('mouse_model_messages.csv','dream_repair_messages.csv')]:
        first, second = pd.read_csv(r/original), pd.read_csv(r/incoming)
        mask = pd.MultiIndex.from_frame(first[['spec','model']]).isin(targets)
        pd.concat([first.loc[~mask],second],ignore_index=True).to_csv(r/original,index=False)


def main():
    """在独立新目录完成复核；0为全部可用，3为仍有失败，2为执行或核验错误。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-work', required=True, type=Path)
    parser.add_argument('--work', required=True, type=Path)
    parser.add_argument('--rscript', default=shutil.which('Rscript') or 'Rscript')
    parser.add_argument('--library', required=True, type=Path)
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    require(platform.system() == 'Linux', '生产入口只支持Linux')
    source, work = args.source_work.resolve(), args.work.resolve()
    require(source != work and not source.is_relative_to(work) and not work.is_relative_to(source), '新旧目录必须分开且互不嵌套')
    require(not work.exists(), '输出目录已存在，请更换新目录；不要覆盖原结果')
    require(args.workers >= 1, 'workers须为正整数')
    verify_release(ROOT)
    import fcntl
    with (source / 'run.lock').open('r') as source_lock:
        fcntl.flock(source_lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
        for folder in ('results','logs','stages','prepared','raw','repair_details'):
            (work/folder).mkdir(parents=True,exist_ok=True)
        try:
            cfg, verified, failures = verify_source(source)
            cfg['analysis_version'] = VERSION; cfg['workers'] = args.workers
            dump(work/'run_config.json', cfg)
            dump(work/'owner.json', dict(signature=fingerprint(cfg),source='targeted_dream_repair'))
            dump(work/'status.json',dict(status='running',task='dream_repair',version=VERSION))
            import_reviewed(source,work,verified)
            targets = failures[['spec','model']].drop_duplicates()
            dump(work/'repair_plan.json',dict(targets=targets.to_dict('records'),
                failures=failures[['spec','model','gene_symbol']].to_dict('records'),
                policy='bobyqa_primary_nloptwrap_check_no_pvalue_selection'))
            print('只复核以下模型：',targets.to_dict('records'),flush=True)
            command([sys.executable,'-m','unittest','discover','-s',ROOT/'tests','-v'],work/'logs/tests.log')
            rrun('preflight.R',['check',args.library.resolve(),ROOT,work/'logs'],work,args.rscript)
            rrun('repair_dream.R',[work,args.library.resolve(),ROOT],work,args.rscript)
            pairs = pd.MultiIndex.from_frame(targets)
            merge_repaired(work,pairs)
            # 重新生成受影响的图表，但模拟、校准和其它模型的数据不重算。
            from pipeline.figures import generate
            generate(work)
            audit = assess(work)
            changed = {n:dict(before=digest,after=sha(work/n)) for n,digest in verified.items()
                       if sha(work/n) != digest}
            allowed = {'results/'+n for n in ['mouse_background_fits.csv.gz','mouse_fit_failures.csv',
                'mouse_model_status.csv','mouse_results.csv.gz','mouse_dream_diagnostics.csv','mouse_model_messages.csv']}
            require(set(changed)<=allowed,'修改超出定点修复范围')
            dump(work/'repair_changes.json',dict(changed=changed,unchanged_files=len(verified)-len(changed),
                    source=str(source),note='方差组分来自原独立提取流程，未冒充本次dream重拟合'))
            dump(work/'status.json',dict(status=audit['status'],task='dream_repair',version=VERSION,audit=audit))
            code = 3 if audit['status'] != 'completed' else 0
        except Exception as exc:
            (work/'logs/error.txt').write_text(traceback.format_exc(),encoding='utf-8')
            dump(work/'status.json',dict(status='failed',task='dream_repair',reason=str(exc)))
            print('停止：',exc,flush=True); code = 2
        finally:
            bundle(work)
        print('结果目录：',work,flush=True)
        return code


if __name__ == '__main__':
    sys.exit(main())
