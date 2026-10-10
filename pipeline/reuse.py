"""复用本流程已经准备好的输入；不下载、不修改来源目录、不复用旧模型结果。"""
import json
import os
import shutil
from pathlib import Path
from .io import dump, require, sha

INPUTS = {
    'panel': ['results/panel.csv', 'results/panel_discovery.csv.gz'],
    'orthologs': ['results/ortholog_panel.csv'],
    'prepare_human': ['prepared/human/metadata.csv', 'prepared/human/expression.csv.gz'],
    'prepare_mouse': ['prepared/mouse_units.csv', 'prepared/mouse_whole_counts.csv.gz',
                      'prepared/mouse_cortex_counts.csv.gz'],
}
INPUT_PARAMETERS = ['seaad_release', 'seaad_expected_donors', 'minimum_nuclei',
                    'minimum_cortex_spots', 'cortex_labels', 'panel_fdr']


def checked_path(root, relative):
    """拒绝绝对路径和越界收据，避免意外引用工作目录外的文件。"""
    require(not Path(relative).is_absolute(), '收据包含绝对路径')
    p = root / relative
    require(p.resolve().is_relative_to(root), f'收据路径越界：{relative}')
    require(p.is_file(), f'输入缺失：{p}')
    return p


def collect(source, config, stage_names=None):
    """依据阶段收据核验必需输入，读取收据内可选的小型元数据。"""
    source = Path(source).resolve()
    owner = json.loads((source / 'owner.json').read_text(encoding='utf-8'))
    require(owner.get('source') == 'official_public_downloads', '来源不是本流程的公开数据工作目录')
    previous = json.loads((source / 'run_config.json').read_text(encoding='utf-8'))
    require(all(previous.get(k) == config.get(k) for k in INPUT_PARAMETERS), '数据准备参数不一致，不能直接复用聚合输入')
    files, receipts = {}, {}
    chosen = list(INPUTS) if stage_names is None else list(stage_names)
    require(set(chosen) <= set(INPUTS), '未知的数据准备阶段')
    for stage in chosen:
        required = INPUTS[stage]
        receipt = source / 'stages' / f'{stage}.json'
        record = json.loads(receipt.read_text(encoding='utf-8'))
        require(record.get('stage') == stage, f'阶段收据身份错误：{stage}')
        outputs = record.get('outputs', {})
        require(set(required) <= set(outputs), f'{stage}收据未登记全部必要输入')
        selected = set(required)
        if stage == 'prepare_human':
            selected.update(p for p in outputs if p.startswith('prepared/human/')
                            and p.endswith(('.csv', '.json')) and '/designs/' not in p)
        for relative in sorted(selected):
            p = checked_path(source, relative)
            print('核验现有输入：', relative, flush=True)
            digest = sha(p)
            require(digest == outputs[relative], f'输入SHA-256不一致：{relative}')
            files[relative] = dict(sha256=digest, bytes=p.stat().st_size)
        receipts[stage] = sha(receipt)
    return dict(source_owner=owner, source_directory=str(source),
                files=files, source_receipts=receipts)


def import_inputs(source, work, config, stage_names=None):
    """先核验再复制，计数优先硬链接以节省空间，模型和图表始终重新生成。"""
    source, work = Path(source).resolve(), Path(work).resolve()
    require(source != work and not work.is_relative_to(source) and not source.is_relative_to(work),
            '来源与新输出必须是互不嵌套的不同目录')
    receipt = work / 'input_reuse.json'
    if receipt.exists():
        prior = json.loads(receipt.read_text(encoding='utf-8'))
        require(prior['source_directory'] == str(source), '续跑不能更换输入来源')
        for relative, info in prior['files'].items():
            p = work / relative
            require(p.is_file() and sha(p) == info['sha256'], f'复用输入已改变：{relative}')
        return prior
    # 来源也上锁，避免读取另一进程正在改写的数据。该函数的生产入口只在Linux使用。
    lock = (source / 'run.lock').open('a')
    try:
        if os.name == 'posix':
            import fcntl
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        record = collect(source, config, stage_names)
        for relative, info in record['files'].items():
            src, dst = source / relative, work / relative
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.exists():
                require(sha(dst) == info['sha256'], f'新目录存在不同内容：{relative}')
                continue
            if relative.endswith('_counts.csv.gz') or relative.endswith('/expression.csv.gz'):
                try:
                    os.link(src, dst)
                    continue
                except OSError:
                    pass
            shutil.copy2(src, dst)
        dump(receipt, record)
        for stage in record['source_receipts']:
            required = INPUTS[stage]
            paths = set(required)
            if stage == 'prepare_human':
                paths.update(p for p in record['files'] if p.startswith('prepared/human/'))
            dump(work / 'stages' / f'{stage}.json', dict(stage=stage,
                 outputs={p: record['files'][p]['sha256'] for p in paths},
                 reused_input=True, source_receipt=record['source_receipts'][stage]))
        print('现有输入核验完成；不启动下载和重新聚合。', flush=True)
        return record
    finally:
        lock.close()
