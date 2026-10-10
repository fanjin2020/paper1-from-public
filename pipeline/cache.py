"""完整入口的数据优先级：合格聚合输入、已下载原始文件、缺失文件下载。"""
import json
import os
from pathlib import Path
import shutil

from .io import dump, require, sha
from .reuse import INPUTS, checked_path, import_inputs


def available_stages(source):
    """只列出文件与阶段记录均齐全的准备阶段；内容正确性随后按摘要核验。"""
    found = []
    for stage, names in INPUTS.items():
        if (source/'stages'/f'{stage}.json').is_file() and all((source/n).is_file() for n in names):
            found.append(stage)
    if 'panel' not in found and 'orthologs' in found:
        found.remove('orthologs')
    return found


def seed_raw_cache(source, work):
    """导入已登记下载文件和可续传残片；完整文件可硬链接，残片必须复制。"""
    source, work = Path(source).resolve(), Path(work).resolve()
    raw = source/'raw'
    records = []
    if not raw.is_dir():
        return records
    from .sources import URLS
    expected = [raw/name for name in URLS]
    expected.extend((raw/'SEAAD').rglob('*.h5ad'))
    for p in expected:
        require(not p.is_file() or p.with_name(p.name+'.download.json').is_file(),
                f'已有文件缺少来源与完整性记录，未重新下载：{p}')
    for receipt in sorted(raw.rglob('*.download.json')):
        relative = receipt.relative_to(source).as_posix().removesuffix('.download.json')
        p = source/relative
        if not p.is_file():
            continue
        p = checked_path(source,relative)
        info = json.loads(receipt.read_text(encoding='utf-8'))
        digest = sha(p)
        require(info.get('bytes') == p.stat().st_size and info.get('sha256') == digest and info.get('url'),
                f'已有原始数据核验失败，未覆盖或重新下载：{relative}')
        dst = work/relative; dst.parent.mkdir(parents=True,exist_ok=True)
        if dst.exists():
            require(sha(dst) == digest, f'新目录内有不同内容：{relative}')
        else:
            try:
                os.link(p,dst)
            except OSError:
                shutil.copy2(p,dst)
        target_receipt = dst.with_name(dst.name+'.download.json')
        if target_receipt.exists():
            require(json.loads(target_receipt.read_text(encoding='utf-8')) == info, '下载记录冲突')
        else:
            shutil.copy2(receipt,target_receipt)
        records.append(dict(path=relative,sha256=digest,bytes=p.stat().st_size))
    for state in sorted(raw.rglob('*.partial.json')):
        base = state.name.removesuffix('.partial.json')
        part = state.with_name(base+'.part')
        rel = state.relative_to(source)
        dst = work/rel.parent/base
        if dst.exists() or not part.is_file():
            continue
        # 只复用明确登记URL和续传校验符的残片，不按文件大小猜下载来源。
        info = json.loads(state.read_text(encoding='utf-8'))
        if not info.get('url') or not info.get('validator'):
            continue
        dst.parent.mkdir(parents=True,exist_ok=True)
        for src in (state,part):
            target = dst.with_name(src.name)
            if not target.exists():
                shutil.copy2(src,target)
    return records


def prepare_existing_data(source, work, config):
    """为all入口选择已有数据；发现损坏时停止，不用重新下载掩盖数据问题。"""
    source, work = Path(source).resolve(), Path(work).resolve()
    require(source.is_dir(), '指定的已有数据目录不存在，请核对--data-work路径')
    require(source != work and not source.is_relative_to(work) and not work.is_relative_to(source),
            '已有数据与新结果目录必须分开且互不嵌套')
    plan_path = work/'data_plan.json'
    if plan_path.exists():
        prior = json.loads(plan_path.read_text(encoding='utf-8'))
        require(prior['source'] == str(source), '续跑不能更换已有数据目录')
        if prior['reused_stages']:
            import_inputs(source,work,config,prior['reused_stages'])
        return set(prior['reused_stages'])
    stages = available_stages(source)
    if stages:
        import_inputs(source,work,config,stages)
    # 完整聚合输入齐全时，不遍历数百GB原始缓存，也不联系远端目录服务。
    records = []
    if set(stages) != set(INPUTS):
        lock = (source/'run.lock').open('a')
        try:
            if os.name == 'posix':
                import fcntl
                fcntl.flock(lock,fcntl.LOCK_SH | fcntl.LOCK_NB)
            records = seed_raw_cache(source,work)
        finally:
            lock.close()
    dump(plan_path,dict(source=str(source),reused_stages=stages,raw_files=records,
                       policy='prepared_then_registered_raw_then_missing_downloads'))
    print('数据检查：复用准备阶段',stages,'；已登记原始文件',len(records),flush=True)
    return set(stages)
