"""下载和运行记录。仅允许使用当前工作目录内登记的文件。"""
import hashlib
import json
import os
import re
from pathlib import Path
import subprocess
import time
import urllib.request
from datetime import datetime, timezone


def require(ok, message):
    """输入不满足分析条件时停止，不改变数据或模型来绕过错误。"""
    if not ok:
        raise ValueError(message)


def sha(path):
    """分块计算内容摘要，不把大文件一次载入内存。"""
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            value.update(block)
    return value.hexdigest()


def dump(path, obj):
    """通过临时文件写入JSON，中断时不留下半个状态文件。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temp.replace(path)


def download(url, path, expected_bytes=None):
    """下载到当前工作目录，保留TLS验证，以ETag或Last-Modified保护断点续传。

    已完成文件必须有本次下载生成的URL、长度和SHA-256记录。远端无摘要时，本地
    SHA-256用于后续完整性检查，不表述为出版方提供的真实性认证。
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    receipt = path.with_name(path.name + '.download.json')
    if path.exists():
        require(receipt.is_file(), f'存在未登记文件，不能作为新下载缓存：{path}')
        record = json.loads(receipt.read_text())
        require(record['url'] == url and record['sha256'] == sha(path)
                and record['bytes'] == path.stat().st_size, f'下载缓存不一致：{path}')
        require(expected_bytes is None or path.stat().st_size == expected_bytes, '远端清单大小改变')
        print(f'已有数据核验通过，跳过下载：{path.name}', flush=True)
        return path
    part = path.with_name(path.name + '.part')
    state = path.with_name(path.name + '.partial.json')
    for attempt in range(1, 7):
        try:
            prior = json.loads(state.read_text()) if state.exists() else {}
            offset = part.stat().st_size if part.exists() else 0
            if offset and (prior.get('url') != url or not prior.get('validator')):
                part.unlink()
                offset = 0
            headers = {'User-Agent': 'paper1-public-reanalysis/1.0', 'Accept-Encoding': 'identity'}
            if offset:
                headers.update(Range=f'bytes={offset}-', **{'If-Range': prior['validator']})
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=180) as response:
                if response.status == 200:
                    offset = 0
                require(response.status in (200, 206), f'HTTP状态{response.status}')
                validator = response.headers.get('ETag') or response.headers.get('Last-Modified')
                if offset:
                    require(response.headers.get('Content-Range', '').startswith(f'bytes {offset}-')
                            and validator == prior['validator'], '远端续传内容或范围改变')
                size = response.headers.get('Content-Length')
                total = int(size) + offset if size else None
                if response.status == 206:
                    match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
                    require(match is not None and int(match[1]) == offset, '无法确认分段响应范围')
                    total = int(match[3])
                dump(state, dict(url=url, validator=validator))
                print(f'下载 {path.name} 第{attempt}次，起点{offset}，总字节{total}', flush=True)
                last = time.monotonic()
                with part.open('ab' if offset else 'wb') as target:
                    while True:
                        block = response.read(1024**2)
                        if not block:
                            break
                        target.write(block)
                        if time.monotonic() - last > 30:
                            print(f'  {path.name}: {target.tell()/1024**3:.2f} GiB', flush=True)
                            last = time.monotonic()
                require(total is None or part.stat().st_size == total, '下载未完成')
                require(expected_bytes is None or part.stat().st_size == expected_bytes, '清单与文件长度不一致')
            record = dict(url=url, bytes=part.stat().st_size, sha256=sha(part), validator=validator,
                          retrieved_utc=datetime.now(timezone.utc).isoformat())
            # 先写收据再移动文件；重启时残片仍可重下，不接受无收据的完整文件。
            dump(receipt, record)
            part.replace(path)
            return path
        except Exception as exc:
            print(f'下载失败 {path.name}: {exc}', flush=True)
            if attempt == 6:
                raise
            time.sleep(min(60, 2**attempt))


def command(args, log, env=None):
    """实时记录外部进程日志，并将非零退出传递给主程序。"""
    log = Path(log)
    log.parent.mkdir(parents=True, exist_ok=True)
    print('RUN', ' '.join(map(str, args)), flush=True)
    with log.open('w', encoding='utf-8') as handle:
        process = subprocess.Popen(list(map(str, args)), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   text=True, encoding='utf-8', errors='replace', env=env or os.environ.copy())
        try:
            for line in process.stdout:
                print(line, end='', flush=True)
                handle.write(line)
                handle.flush()
            rc = process.wait()
        except BaseException:
            process.terminate()
            process.wait()
            raise
    require(rc == 0, f'进程退出码{rc}；查看{log}')
