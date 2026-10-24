"""Import only integrity-checked local copies of public download objects."""
import contextlib
import json
import os
from pathlib import Path
import re
import shutil
from urllib.parse import urlparse

from .io import dump, require, sha


def destination(url, release):
    from .sources import URLS
    exact = {value: name for name, value in URLS.items()}
    if url in exact:
        return Path(exact[url])
    parsed = urlparse(url)
    suffix = "_SEAAD_DFC_RNAseq_final-nuclei." + release + ".h5ad"
    if (parsed.scheme == "https" and parsed.netloc == "sea-ad-single-cell-profiling.s3.amazonaws.com"
            and parsed.path.startswith("/PFC/RNAseq/donors_objects/")
            and not parsed.query and not parsed.fragment and parsed.path.endswith(suffix)):
        name = Path(parsed.path).name
        require(re.fullmatch(r"[A-Za-z0-9_.-]+", name) is not None, "Unsafe download name")
        return Path("SEAAD") / name
    return None


@contextlib.contextmanager
def read_lock(root):
    lock = root / "run.lock"
    if not lock.exists() or os.name != "posix":
        yield
        return
    import fcntl
    with lock.open("rb") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
        yield


def seed_raw_cache(source, work, config):
    """Never import prepared matrices, NPZ, stage receipts, fits or simulations."""
    source, work = Path(source).resolve(), Path(work).resolve()
    require(source.is_dir(), "Raw cache directory not found: " + str(source))
    require(source != work and not source.is_relative_to(work) and not work.is_relative_to(source),
            "Raw cache and new core output must be separate, non-nested directories.")
    raw = source / "raw" if (source / "raw").is_dir() else source
    records = []
    with read_lock(source):
        for receipt in sorted(raw.rglob("*.download.json")):
            info = json.loads(receipt.read_text(encoding="utf-8-sig"))
            relative = destination(info.get("url", ""), config["seaad_release"])
            if relative is None:
                continue
            path = receipt.with_name(receipt.name.removesuffix(".download.json"))
            if not path.is_file():
                continue
            require(path.resolve().is_relative_to(raw.resolve()), "Download path escapes cache")
            print("VERIFY RAW", path, flush=True)
            digest = sha(path)
            require(info.get("sha256") == digest and info.get("bytes") == path.stat().st_size,
                    "Raw cache is damaged; left untouched: " + str(path))
            target = work / "raw" / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target_receipt = target.with_name(target.name + ".download.json")
            if target.exists():
                require(sha(target) == digest, "Conflicting raw versions: " + str(target))
            else:
                try:
                    os.link(path, target)
                except OSError:
                    shutil.copy2(path, target)
            dump(target_receipt, info)
            records.append({"source": str(path), "destination": relative.as_posix(),
                            "url": info["url"], "bytes": info["bytes"], "sha256": digest})
    return records


def import_raw_sources(sources, work, config):
    work = Path(work)
    report = work / "raw_cache_import.json"
    canonical = [str(Path(p).expanduser().resolve()) for p in sources]
    if report.exists():
        previous = json.loads(report.read_text(encoding="utf-8"))
        require(previous.get("sources") == canonical, "Same-run continuation cannot change raw-cache sources.")
        return previous["files"]
    records = []
    for source in canonical:
        records.extend(seed_raw_cache(source, work, config))
    dump(report, {"policy": "download_objects_only", "sources": canonical, "files": records,
                  "prepared_inputs_imported": 0, "previous_results_imported": 0})
    return records
