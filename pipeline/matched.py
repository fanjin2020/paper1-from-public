#!/usr/bin/env python3
"""Matched human sensitivity, seven mouse model arms and count-to-inference validation."""

import argparse
import concurrent.futures as cf
import contextlib
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import statistics
import subprocess
import sys
import threading
import time
import traceback

from pipeline.configuration import VERSION
ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "analysis" / "matched"
REQUIRED = {
    "panel": ["results/panel.csv", "results/panel_discovery.csv.gz"],
    "orthologs": ["results/ortholog_panel.csv"],
    "prepare_human": ["prepared/human/metadata.csv", "prepared/human/expression.csv.gz"],
    "prepare_mouse": ["prepared/mouse_units.csv", "prepared/mouse_whole_counts.csv.gz",
                      "prepared/mouse_cortex_counts.csv.gz"],
}
CHILDREN = set()
GUARD = threading.Lock()
STOP = threading.Event()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8-sig"))
    except (OSError, ValueError):
        return {}


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    os.replace(tmp, path)


@contextlib.contextmanager
def lock_output(out):
    p = Path(out) / "run.lock"
    with p.open("a+b") as f:
        try:
            if os.name == "posix":
                import fcntl
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
            else:
                import msvcrt
                f.seek(0)
                if not f.read(1):
                    f.write(b"0")
                    f.flush()
                f.seek(0)
                msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 1)
        except OSError as e:
            raise RuntimeError("Output is locked by another run; do not delete run.lock.") from e
        yield


def resources():
    cpu = os.cpu_count() or 1
    if hasattr(os, "sched_getaffinity"):
        cpu = len(os.sched_getaffinity(0))
    available = None
    if Path("/proc/meminfo").exists():
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                available = int(line.split()[1]) * 1024
    load = os.getloadavg()[0] if hasattr(os, "getloadavg") else 0.
    return {"cpus": cpu, "available_GiB": available / 2**30 if available else None,
            "load_1min": load}


@contextlib.contextmanager
def source_read_lock(source):
    # Never create or delete a lock in the original analysis directory.
    path = source / "run.lock"
    if os.name != "posix" or not path.exists():
        yield
        return
    import fcntl
    with path.open("rb") as f:
        try:
            fcntl.flock(f, fcntl.LOCK_SH | fcntl.LOCK_NB)
        except OSError as e:
            raise RuntimeError("The original data work is still locked by a writer; let it finish first.") from e
        yield


def worker_limit(requested):
    r = resources()
    if requested:
        require(requested <= r['cpus'], 'Requested workers exceed the available CPU affinity')
        selected = requested
        if r['available_GiB'] is not None and r['available_GiB'] < 24 + 2 * selected:
            print('WARNING: limited available memory; monitor memory while running', flush=True)
    else:
        memory = max(1, int((r['available_GiB'] - 24) / 4)) if r['available_GiB'] else 36
        selected = max(1, min(36, memory, r['cpus'] - math.ceil(r['load_1min']) - 4))
    print('RESOURCES', json.dumps(r), 'WORKERS', selected, flush=True)
    return selected


def discover_source(args, out):
    require(args.data_work is not None, 'Specify --data-work with prepared core inputs')
    source = args.data_work.expanduser().resolve()
    config = read_json(source / 'run_config.json')
    require(config.get('minimum_nuclei', 0) > 0, 'Missing input preparation settings')
    records = source_fingerprint(source)
    for name, expected in records.items():
        require(sha(checked_relative(source, name)) == expected, 'Input checksum differs: ' + name)
    require(config.get('ebayes_robust') is False and config.get('dream_ebayes_robust') is False,
            'Primary models require ordinary eBayes; robust moderation is a separate sensitivity arm')
    inputs = {'prepared': records, 'configuration': {'run_config.json': sha(source / 'run_config.json')}}
    scientific = {k: config[k] for k in ('minimum_nuclei','minimum_cortex_spots',
                  'cortex_labels','panel_fdr','ebayes_robust','dream_ebayes_robust')}
    return source, inputs, scientific


def checked_relative(root, relative):
    p = (Path(root) / relative).resolve()
    require(not Path(relative).is_absolute() and p.is_relative_to(Path(root).resolve()),
            "Receipt path escapes source: " + relative)
    require(p.is_file(), "Missing input: " + str(p))
    return p


def source_fingerprint(root):
    records = {}
    for stage, names in REQUIRED.items():
        receipt = read_json(root / "stages" / (stage + ".json"))
        require(receipt.get("stage") == stage, "Missing source receipt: " + stage)
        for name in names:
            require(name in receipt.get("outputs", {}), "Input missing from source receipt: " + name)
            records[name] = receipt["outputs"][name]
    return records




def r_environment(args, out):
    exe = args.rscript or shutil.which("Rscript")
    require(exe and Path(exe).is_file(), "Rscript not found; use --rscript /path/to/Rscript")
    # Bypass the Windows architecture launcher so cancellation targets the R process itself.
    if os.name == "nt" and (Path(exe).parent / "x64" / "Rscript.exe").is_file():
        exe = str(Path(exe).parent / "x64" / "Rscript.exe")
    local = out / "runtime" / "R_library"
    local.mkdir(parents=True, exist_ok=True)
    libraries = [local, *[Path(p).expanduser().resolve() for p in args.r_library]]
    libraries = list(dict.fromkeys(str(p) for p in libraries if p.is_dir()))
    env = os.environ.copy()
    env["PAPER1_R_LIBRARIES"] = os.pathsep.join(libraries)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    env["BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS"] = "FALSE"
    log = out / "logs" / "R_check.log"
    run_command([exe, "--vanilla", str(HERE / "analysis.R"), "check",
                 str(out / "runtime" / "environment.json"), str(local),
                 "TRUE" if args.setup and not args.offline else "FALSE"], env, log)
    return str(Path(exe).resolve()), env, read_json(out / "runtime" / "environment.json")


def stop_children():
    STOP.set()
    with GUARD:
        procs = list(CHILDREN)
    for p in procs:
        if p.poll() is None:
            try:
                if os.name == "posix":
                    os.killpg(p.pid, signal.SIGTERM)
                else:
                    p.terminate()
            except ProcessLookupError:
                pass
    for p in procs:
        try:
            p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            if os.name == "posix":
                os.killpg(p.pid, signal.SIGKILL)
            else:
                p.kill()
            p.wait()


def run_command(command, env, log):
    log.parent.mkdir(parents=True, exist_ok=True)
    require(not STOP.is_set(), "Interrupted")
    with log.open("w", encoding="utf-8") as f:
        p = subprocess.Popen(list(map(str, command)), stdout=f, stderr=subprocess.STDOUT, env=env,
                             start_new_session=os.name == "posix")
        with GUARD:
            CHILDREN.add(p)
        try:
            code = p.wait()
        except BaseException:
            stop_children()
            raise
        finally:
            with GUARD:
                CHILDREN.discard(p)
    if code:
        print(log.read_text(encoding="utf-8", errors="replace")[-6500:], flush=True)
        raise RuntimeError("R task failed; read " + str(log))


def task_health(folder, job):
    """File integrity alone cannot certify an entirely failed analysis arm."""
    import csv
    import gzip
    phase = job.get("phase")
    if phase == "mouse":
        with gzip.open(folder / "background_fits.csv.gz", "rt", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        expected = {"baseline", "cage_fixed", "cage_correlation", "dream_shared", "dream_native",
                    "mundlak_shared", "mundlak_native"}
        require({r["model"] for r in rows} == expected, "Missing empirical mouse model: " + job["id"])
        available = {m: 0 for m in expected}
        identities = set()
        for row in rows:
            key = (row["model"], row["gene_symbol"])
            require(key not in identities, "Duplicate mouse fit: " + str(key))
            identities.add(key)
            if row["status"] == "available":
                require(math.isfinite(float(row["p"])), "Nonfinite available mouse fit")
                available[row["model"]] += 1
        require(all(available.values()), "Whole empirical mouse model failed: " + str(available))
    elif phase == "counts":
        with (folder / "metrics.csv").open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        require(len(rows) == 5 * (job["end"] - job["start"]), "Incomplete count task: " + job["id"])
        keys = [(r["scenario"], r["replicate"], r["model"]) for r in rows]
        require(len(set(keys)) == len(keys), "Duplicate count dataset/model")
        require(all(r["usable_run"].lower() in ("true", "1") for r in rows),
                "Whole count model failed: " + job["id"])
    diagnostics = folder / "diagnostics.csv"
    if diagnostics.is_file():
        with diagnostics.open(encoding="utf-8") as stream:
            require(not any(r.get("status") == "retry_merge_failed" for r in csv.DictReader(stream)),
                    "Optimizer retry merge failure: " + job["id"])


def summary_ok(out):
    summary = read_json(out / "results" / "analysis_summary.json")
    require(summary.get("status") == "ANALYSES_COMPLETE",
            "Analysis summary requires review: " + str(summary.get("status", "missing")))
    return summary


def receipt_ok(out, job, protocol_hash, verify_hash=True, validate_health=True):
    r = read_json(out / "receipts" / (job["id"] + ".json"))
    if not r:
        return False
    require(r.get("protocol_hash") == protocol_hash and r.get("job") == job,
            "Receipt belongs to different protocol: " + job["id"])
    for name, v in r.get("files", {}).items():
        p = checked_relative(out, name)
        require(p.is_file() and p.stat().st_size == v["bytes"] and
                (not verify_hash or sha(p) == v["sha256"]),
                "Damaged completed task: " + job["id"] + "; preserve output and inspect before resuming")
    require(bool(r.get("files")), "Empty completion receipt")
    if validate_health:
        task_health(out / "tasks" / job["id"], job)
    return True


def run_job(job, out, protocol, exe, env, cfg):
    if STOP.is_set():
        return None
    if receipt_ok(out, job, protocol["hash"]):
        return "reused"
    start = time.monotonic()
    folder = out / "tasks" / job["id"]
    folder.mkdir(parents=True, exist_ok=True)
    request = dict(job, output=str(folder), root=str(out), config=str(cfg))
    atomic_json(out / "jobs" / (job["id"] + ".json"), request)
    run_command([exe, "--vanilla", HERE / "analysis.R", "job", out / "jobs" / (job["id"] + ".json")],
                env, out / "logs" / (job["id"] + ".log"))
    done = read_json(folder / "complete.json")
    require(done.get("id") == job["id"], "Task did not publish completion: " + job["id"])
    paths = [folder / n for n in done.get("files", [])]
    require(paths and all(p.is_file() and p.resolve().is_relative_to(folder.resolve()) for p in paths),
            "Incomplete task files: " + job["id"])
    task_health(folder, job)
    files = {p.relative_to(out).as_posix(): {"sha256": sha(p), "bytes": p.stat().st_size} for p in paths}
    atomic_json(out / "receipts" / (job["id"] + ".json"),
                {"job": job, "protocol_hash": protocol["hash"], "seconds": time.monotonic() - start,
                 "finished": time.time(), "files": files})
    return "completed"


def make_jobs(meta, args):
    jobs = []
    for i, name in enumerate(meta["contexts"]):
        jobs.append({"id": f"human_{i:03d}", "phase": "human", "context": name, "index": i + 1})
    for spec in meta["mouse_specs"]:
        jobs.append({"id": "mouse_" + spec, "phase": "mouse", "spec": spec})
    for s in range(8):
        for lo in range(0, args.reps, args.batch_size):
            jobs.append({"id": f"counts_s{s:02d}_r{lo:05d}", "phase": "counts", "scenario": s,
                         "start": lo, "end": min(args.reps, lo + args.batch_size)})
    return jobs






def run_phase(jobs, concurrency, out, protocol, exe, env, cfg):
    pending = [j for j in jobs if not receipt_ok(out, j, protocol["hash"])]
    print("PHASE", jobs[0]["phase"] if jobs else "empty", "REMAINING", len(pending), "WORKERS", concurrency, flush=True)
    executor = cf.ThreadPoolExecutor(max_workers=concurrency)
    futures = {executor.submit(run_job, j, out, protocol, exe, env, cfg): j for j in pending}
    try:
        for i, f in enumerate(cf.as_completed(futures), 1):
            f.result()
            print("DONE", futures[f]["id"], f"{i}/{len(pending)}", flush=True)
    except BaseException:
        stop_children()
        for f in futures:
            f.cancel()
        raise
    finally:
        executor.shutdown(wait=True, cancel_futures=True)


def status(out, full=False):
    p = read_json(out / "protocol.json")
    print("STATE", read_json(out / "state.json").get("status", "NOT_STARTED"))
    lock = out / "run.lock"
    if os.name == "posix" and lock.is_file():
        import fcntl
        with lock.open("rb") as f:
            try:
                fcntl.flock(f, fcntl.LOCK_SH | fcntl.LOCK_NB)
                print("WRITER LOCK IDLE")
            except BlockingIOError:
                print("WRITER LOCK ACTIVE")
    if not p:
        print("No protocol yet: dependency check or input preparation may still be running.")
        return
    for phase in ("human", "mouse", "counts"):
        jobs = [j for j in p["jobs"] if j["phase"] == phase]
        finished, seconds, repeats = 0, [], 0
        for job in jobs:
            if receipt_ok(out, job, p["hash"], verify_hash=full):
                finished += 1
                receipt = read_json(out / "receipts" / (job["id"] + ".json"))
                seconds.append(receipt["seconds"])
                repeats += job.get("end", 1) - job.get("start", 0)
        print(f"{phase}: {finished}/{len(jobs)} tasks ({100*finished/max(1,len(jobs)):.1f}%), completed datasets/units={repeats}")
        if (phase == "counts" and len(seconds) >= 5 and finished < len(jobs)
                and read_json(out / "state.json").get("status") == "RUNNING"):
            w = read_json(out / "state.json").get("workers", 1)
            estimate = statistics.median(seconds) * (len(jobs) - finished) / w / 3600
            print(f"Approximate remaining count-simulation time: {estimate:.1f} h; excludes other phases and summaries.")
    print("Receipt checks:", "SHA256" if full else "identity and sizes; add --verify-status for SHA256")




def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, default=Path.cwd())
    p.add_argument("--data-work", type=Path)
    p.add_argument("--output", type=Path, default=Path("results/matched"))
    p.add_argument("--rscript")
    p.add_argument("--r-library", action="append", default=[])
    p.add_argument("--workers", type=int, default=0, help="0=resource-aware up to 36; explicit count is honored within CPU affinity")
    p.add_argument("--reps", type=int, default=200, help="Independent datasets per count scenario; default 8 x 200")
    p.add_argument("--batch-size", type=int, default=2)
    p.add_argument("--sim-genes", type=int, default=0, help="0=all common expressed mouse genes; nonzero is reduced-scope validation")
    p.add_argument("--seed", type=int, default=20261015)
    p.add_argument("--setup", action="store_true")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--smoke", action="store_true", help="Small implementation test, never manuscript evidence")
    p.add_argument("--status", action="store_true")
    p.add_argument("--verify-status", action="store_true")
    p.add_argument("--resume", action="store_true", help="Continue this same raw-input experiment only")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    args.output = args.output.expanduser().resolve()
    args.project_root = args.project_root.expanduser().resolve()
    out = args.output
    if args.status or args.verify_status:
        status(out, args.verify_status)
        return 0
    if args.smoke:
        args.reps, args.batch_size, args.sim_genes = 2, 1, 300
    require(args.reps >= 2 and args.batch_size >= 1 and args.workers >= 0, "Invalid repetitions or workers")
    require(0 <= args.seed and args.seed + 100003 * (args.reps - 1) + 1009 * 7 < 2**31,
            "Seed/repetition combination exceeds the explicit R integer seed range")
    require(args.sim_genes == 0 or args.sim_genes >= 100, "Need at least 100 genes for moderation")
    require(args.data_work is not None, "Use run.py to supply this run's freshly prepared core")
    require(not out.exists() or not any(out.iterdir()) or
            read_json(out / "owner.json").get("runner") == "paper1_matched",
            "Output is nonempty and not owned by this runner")
    out.mkdir(parents=True, exist_ok=True)
    with lock_output(out), contextlib.ExitStack() as source_locks:
        from pipeline.lifecycle import claim, current_source, parent_run, code_identity
        require(parent_run() is not None and args.data_work.resolve() == out.parent / 'core',
                "Use run.py: downstream inputs must belong to this composite raw-input experiment.")
        source_id = current_source(args.data_work, parent_run())
        run_id = claim(out, 'matched', args.resume, source_id, code_identity(ROOT))
        atomic_json(out / "owner.json", {"runner": "paper1_matched", "version": VERSION})
        try:
            workers = worker_limit(args.workers)
            atomic_json(out / "state.json", {"status": "PREPARING", "pid": os.getpid(), "workers": workers})
            source, inputs, source_config = discover_source(args, out)
            require(source != out and not source.is_relative_to(out) and not out.is_relative_to(source),
                    "Source and output must not be nested")
            source_locks.enter_context(source_read_lock(source))
            exe, env, versions = r_environment(args, out)
            code = {p.name: sha(p) for p in HERE.glob("*.R")}
            code.update({"../../pipeline/matched.py": sha(Path(__file__)),
                         "../common.R": sha(HERE.parent / "common.R"),
                         "../diagnostics.R": sha(HERE.parent / "diagnostics.R")})
            settings = {"version": VERSION, "run_id": run_id, "code": code, "inputs": inputs["prepared"], "R": versions,
                        "reps": args.reps, "batch_size": args.batch_size, "seed": args.seed,
                        "sim_genes": args.sim_genes, "smoke": args.smoke,
                        "source_config": source_config,
                        "empirical_mundlak": "arithmetic cage means of all baseline columns; two weight policies; no weighted identity assumed",
                        "simulation": "NB gene-dispersion calibration; gene-specific latent log2 cage SD and genotype-related cage mean; re-filter/TMM/voom/eBayes/BH",
                        "simulation_scope": "whole-section common-background genes; same 74-animal allocation; 8 prespecified scenarios"}
            digest = signature(settings)
            previous = read_json(out / "protocol.json")
            require(not previous or previous["hash"] == digest,
                    "Scientific settings, inputs, R environment or source changed. Use a new output; do not erase old receipts.")
            atomic_json(out / "input_manifest.json", {"source": str(source), **inputs})
            cfg = out / "analysis_config.json"
            # Parallelize independent designs, not nested SOCK backends sharing a port.
            inner = 1
            atomic_json(cfg, {"source": str(source), "out": str(out), "minimum_nuclei": source_config["minimum_nuclei"],
                              "smoke": args.smoke, "sim_genes": args.sim_genes, "seed": args.seed,
                              "reps": args.reps, "inner_workers": inner, "hash": digest})
            prep = read_json(out / "runtime" / "prepared_receipt.json")
            if not prep or prep.get("hash") != digest:
                run_command([exe, "--vanilla", HERE / "analysis.R", "prepare", cfg], env, out / "logs" / "prepare.log")
                paths = sorted((out / "cache").glob("*.rds"))
                require(paths, "Preparation produced no analysis caches")
                paths += sorted((out / "results").glob("*.csv"))
                paths += [out / "runtime" / "prepared_manifest.json"]
                atomic_json(out / "runtime" / "prepared_receipt.json", {"hash": digest,
                            "files": {p.relative_to(out).as_posix(): sha(p) for p in paths}})
            else:
                for name, expected in prep["files"].items():
                    require((out / name).is_file() and sha(out / name) == expected, "Prepared input damaged: " + name)
            meta = read_json(out / "runtime" / "prepared_manifest.json")
            jobs = make_jobs(meta, args)
            protocol = {"hash": digest, "scientific": settings, "jobs": jobs}
            atomic_json(out / "protocol.json", protocol)
            atomic_json(out / "state.json", {"status": "RUNNING", "pid": os.getpid(), "workers": workers,
                                             "started": time.time(), "validation_only": args.smoke})
            for phase, capacity in (("human", min(12, workers)), ("mouse", min(6, workers)), ("counts", workers)):
                run_phase([j for j in jobs if j["phase"] == phase], capacity, out, protocol, exe, env, cfg)
            require(all(receipt_ok(out, j, digest) for j in jobs), "Not all tasks completed")
            run_command([exe, "--vanilla", HERE / "analysis.R", "summarise", cfg], env, out / "logs" / "summarise.log")
            summary_ok(out)
            for name, expected in {**inputs["prepared"], **inputs["configuration"]}.items():
                require(sha(source / name) == expected, "Original input changed during execution")
            for name, expected in code.items():
                require(sha(HERE / name) == expected, "Analysis code changed during execution; do not edit a running source snapshot")
            final_files = [p for folder in ("results", "figures") for p in (out / folder).rglob("*") if p.is_file()]
            atomic_json(out / "summary_receipt.json", {"protocol_hash": digest,
                        "files": {p.relative_to(out).as_posix(): sha(p) for p in final_files}})
            atomic_json(out / "state.json", {"status": "COMPLETED", "workers": workers, "finished": time.time(),
                                             "validation_only": args.smoke, "protocol_hash": digest})
            print("COMPLETED", out, flush=True)
            return 0
        except BaseException as e:
            stop_children()
            atomic_json(out / "state.json", {"status": "INTERRUPTED" if isinstance(e, KeyboardInterrupt) else "FAILED",
                                             "reason": str(e), "time": time.time(), "workers": locals().get("workers")})
            if isinstance(e, KeyboardInterrupt):
                print("Interrupted. Add --resume to continue this same run using verified task receipts.")
                return 130
            raise


if __name__ == "__main__":
    raise SystemExit(main())
