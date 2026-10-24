#!/usr/bin/env python3
"""Paired coefficient calibration and five-model grouped-data simulations.

Run --help for setup, resumable execution, status, and smoke-test options.
Official R model fitting is defined in analysis/benchmark.R.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import contextlib
import csv
import hashlib
import io
import itertools
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import traceback

for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_key] = "1"

from pipeline.configuration import VERSION
PAIRED_SEED = 20261015
CAGE_SEED = 71004000
METHODS = ["baseline", "cage_fixed", "cage_correlation", "dream", "mundlak_dream"]
PY_REQUIREMENTS = ["numpy==2.2.6", "scipy==1.15.3", "pandas==2.2.3", "matplotlib==3.10.1"]

ACTUAL_INPUT = None


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False), encoding="utf-8")
    os.replace(tmp, path)


def atomic_csv(frame, path):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    frame.to_csv(tmp, index=False, compression="gzip" if path.suffix == ".gz" else None)
    os.replace(tmp, path)


def read_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


@contextlib.contextmanager
def output_lock(out, probe=False):
    path = Path(out) / "run.lock"
    if probe and not path.exists():
        yield True
        return
    if not probe:
        path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("rb" if probe else "a+b") as handle:
        if path.stat().st_size == 0:
            if not probe:
                handle.write(b"0")
                handle.flush()
        handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBRLCK if probe else msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), (fcntl.LOCK_SH if probe else fcntl.LOCK_EX) | fcntl.LOCK_NB)
        except OSError as exc:
            if probe:
                yield False
                return
            raise RuntimeError("Output is locked by another run. Stop that run first; do not delete run.lock.") from exc
        try:
            yield True
        finally:
            handle.seek(0)
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def load_science():
    global np, pd, stats
    import numpy as np
    import pandas as pd
    from scipy import stats


def bootstrap(args):
    missing = []
    for name in ("numpy", "scipy", "pandas", "matplotlib"):
        try:
            __import__(name)
        except ImportError:
            missing.append(name)
    if not missing:
        return
    require(args.setup and not args.offline,
            "Missing Python packages: " + ", ".join(missing) + ". Use --setup (without --offline).")
    import venv
    folder = Path(__file__).resolve().parents[1] / ".benchmark-python"
    exe = folder / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    if not exe.exists():
        venv.EnvBuilder(with_pip=True).create(folder)
    subprocess.run([str(exe), "-m", "pip", "install", *PY_REQUIREMENTS], check=True)
    cmd = [str(exe), str(Path(__file__).resolve()), *sys.argv[1:]]
    raise SystemExit(subprocess.call(cmd))


def basis(x):
    u, s, _ = np.linalg.svd(np.asarray(x, float), full_matrices=False)
    tol = max(x.shape) * np.finfo(float).eps * s[0]
    return u[:, s > tol]


def zscore(x):
    x = np.asarray(x, float)
    require(x.std() > 0, "Constant design column")
    return (x - x.mean()) / x.std()


def geometry(x0, x1, w):
    x0, x1, w = np.asarray(x0, float), np.asarray(x1, float), np.asarray(w, float)
    require(np.all(w > 0) and np.all(np.isfinite(w)), "Invalid weights")
    require(np.linalg.matrix_rank(x0) == x0.shape[1] and
            np.linalg.matrix_rank(x1) == x1.shape[1], "Rank-deficient paired design")
    require(np.linalg.matrix_rank(np.column_stack([x1, x0])) == x1.shape[1], "Models must be nested")
    a0 = np.linalg.solve(x0.T @ (w[:, None] * x0), x0.T * w)
    a1 = np.linalg.solve(x1.T @ (w[:, None] * x1), x1.T * w)
    d = a1[1] - a0[1]
    h = np.sum(x1 * a1.T, axis=1)
    require(len(w) > x1.shape[1] and h.max() < 1 - 1e-10, "No residual information")
    return dict(x0=x0, x1=x1, w=w, a0=a0, a1=a1, d=d, h=h,
                unit=float(np.sum(d * d / w)), df=len(w) - x1.shape[1],
                independent_unit=float(np.sum((a0[1] ** 2 + a1[1] ** 2) / w)))


def generic_design(spec):
    rng = np.random.default_rng(np.random.SeedSequence([PAIRED_SEED, 1000, spec["design_index"]]))
    n = spec["n"]
    age = zscore(rng.normal(size=n))
    base = np.column_stack([np.ones(n), age])
    x = rng.normal(size=n)
    if spec["leverage"] == "high":
        x[0] += 2 * np.sqrt(n)
    x = zscore(x - base @ np.linalg.lstsq(base, x, rcond=None)[0])
    orth = rng.normal(size=(n, 3))
    common = np.column_stack([base, x])
    orth -= common @ np.linalg.lstsq(common, orth, rcond=None)[0]
    q = np.linalg.qr(orth)[0] * np.sqrt(n)
    r = spec["partial_correlation"]
    z = np.column_stack([r * x + np.sqrt(1 - r * r) * q[:, 0], q[:, 1:]])
    x0 = np.column_stack([np.ones(n), x, age])
    w = np.ones(n) if spec["weights"] == "uniform" else np.exp(.4 * np.clip(age, -2, 2))
    w /= w.mean()
    return geometry(x0, np.column_stack([x0, z]), w)


CANCELLATION_RULE = "lowest_index_within_256eps_scaled_minimum_v1"


def cancellation_direction(c):
    """Resolve numerically tied coordinates before projecting into the null space."""
    c = np.asarray(c, float)
    cn = np.linalg.norm(c)
    require(np.isfinite(cn) and cn > 1e-10, "degenerate_difference")
    tolerance = 256 * np.finfo(float).eps * max(1., cn)
    index = int(np.flatnonzero(np.abs(c) <= np.min(np.abs(c)) + tolerance)[0])
    v = np.eye(len(c))[index]
    v -= c * (v @ c) / (cn * cn)
    require(np.linalg.norm(v) > 1e-10, "Degenerate cancellation direction")
    return index, v, tolerance


def means(g):
    z = g["x1"][:, -3:]
    c = g["d"] @ z
    cn = np.linalg.norm(c)
    require(cn > 1e-10 and g["unit"] > 1e-14, "degenerate_difference")
    _, v, _ = cancellation_direction(c)
    cancellation = z @ v
    cancellation *= .5 / np.sqrt(np.mean(cancellation ** 2))
    common = .2 * g["x0"][:, 1] + .1 * g["x0"][:, 2]
    alternative = z @ c * np.sqrt(g["unit"]) / (cn * cn)
    return [("quality_zero", common, 0.),
            ("quality_nonzero_delta_zero", common + cancellation, 0.),
            ("delta_one_reference_SD", common + alternative, np.sqrt(g["unit"])),
            ("delta_two_reference_SD", common + 2 * alternative, 2 * np.sqrt(g["unit"]))]


def wilson(k, n):
    if not n:
        return float("nan"), float("nan")
    z = 1.959963984540054
    p = k / n
    center = (p + z*z/(2*n)) / (1+z*z/n)
    half = z * np.sqrt(p*(1-p)/n+z*z/(4*n*n)) / (1+z*z/n)
    return max(0., center-half), min(1., center+half)


def sufficient(est, se, truth, crit, planned, joint):
    valid = np.isfinite(est) & np.isfinite(se) & (se > 1e-12)
    e, s = est[valid] - truth, se[valid]
    return dict(planned=int(planned), available=int(valid.sum()),
                not_estimable=0, numerical_failed=int(planned-valid.sum()),
                sum_error=float(e.sum()), sum_error2=float(e @ e), sum_SE=float(s.sum()),
                sum_SE2=float(s @ s), covered=int(np.sum(np.abs(e) <= crit*s)),
                rejected=int(np.sum(np.abs(est[valid]) > crit*s)),
                sum_width=float(np.sum(2*crit*s)), joint_F_rejected=int(joint))


def paired_batch(job):
    load_science()
    spec = job["spec"]
    g = generic_design(spec)
    start, end = job["start"], job["end"]
    count = end-start
    fields = {k: v for k, v in spec.items() if k != "design_index"}
    fields.update(n=spec["n"], df=g["df"], max_leverage=float(g["h"].max()),
                  unit_variance=g["unit"], batch_start=start)
    records = []
    if g["unit"] < 1e-14:
        for mean_name in ("quality_zero", "quality_nonzero_delta_zero", "delta_one_reference_SD", "delta_two_reference_SD"):
            for method in ("model_based", "HC3", "independence_comparator"):
                zero = sufficient(np.array([]), np.array([]), 0., 1., count, 0)
                zero.update(numerical_failed=0, not_estimable=count)
                records.append(dict(**fields, mean_case=mean_name, method=method,
                                    true_delta=float("nan") if mean_name.startswith("delta_") else 0.,
                                    oracle_SD=0., status="degenerate_difference", **zero))
        return pd.DataFrame(records)
    error_index = ("matched_gaussian", "heteroskedastic_gaussian", "matched_student_t5").index(spec["error_law"])
    # Per-repetition streams make draws independent of worker count and batches.
    noise = []
    variance = 1/g["w"]
    if error_index == 1:
        mult = np.exp(.7*np.clip(zscore(g["x1"][:, -3]), -2, 2))
        variance = variance * mult/mult.mean()
    for r in range(start, end):
        rng = np.random.default_rng(np.random.SeedSequence([PAIRED_SEED, 2000, spec["design_index"], error_index, r]))
        e = rng.standard_t(5, len(variance))*np.sqrt(3/5) if error_index == 2 else rng.standard_normal(len(variance))
        noise.append(e*np.sqrt(variance))
    noise = np.asarray(noise).T
    oracle_sd = float(np.sqrt(np.sum(g["d"]**2*variance)))
    crit = stats.t.ppf(.975, g["df"])
    for mean_name, mu, truth in means(g):
        require(abs(g["d"] @ mu-truth) < 1e-9, "Constructed coefficient difference failed identity")
        y = mu[:, None]+noise
        e1 = y-g["x1"]@(g["a1"]@y)
        e0 = y-g["x0"]@(g["a0"]@y)
        rss = np.sum(g["w"][:, None]*e1**2, axis=0)
        delta = g["d"]@y
        f = np.maximum(0, (np.sum(g["w"][:, None]*e0**2, axis=0)-rss)/3)/(rss/g["df"])
        joint = int(np.sum(stats.f.sf(f, 3, g["df"]) < .05))
        ses = [("model_based", np.sqrt(rss/g["df"]*g["unit"])),
               ("HC3", np.sqrt((g["d"]**2/(1-g["h"])**2)@(e1**2))),
               ("independence_comparator", np.sqrt(rss/g["df"]*g["independent_unit"]))]
        for method, se in ses:
            records.append(dict(**fields, mean_case=mean_name, method=method, true_delta=float(truth),
                                oracle_SD=oracle_sd, status="regular", **sufficient(delta, se, truth, crit, count, joint)))
    return pd.DataFrame(records)


def synthetic(size, mixed, fraction, seed):
    n = 72
    rng = np.random.default_rng(seed)
    k, m, target = n//size, round(n//size*mixed), round(n*fraction)
    options = [p for p in range(k-m+1) if m <= target-p*size <= m*(size-1)]
    require(options, "Infeasible synthetic allocation")
    pure = options[len(options)//2]
    nums, left = [1]*m, target-pure*size-m
    while left:
        choices = [j for j, v in enumerate(nums) if v < size-1]
        low = min(nums[j] for j in choices)
        j = int(rng.choice([j for j in choices if nums[j] == low]))
        nums[j] += 1
        left -= 1
    nums += [size]*pure+[0]*(k-m-pure)
    rng.shuffle(nums)
    rows = []
    for j, num in enumerate(nums):
        gs = [1]*num+[0]*(size-num)
        rng.shuffle(gs)
        offset = len(rows)
        rows.extend(dict(sample=f"s{offset+i}", cage=j, genotype=v) for i, v in enumerate(gs))
    return pd.DataFrame(rows)


def observed_layout(path):
    data = pd.read_csv(path, dtype={'SAMPLE': str, 'cage': str}).sort_values('SAMPLE').reset_index(drop=True)
    require({'SAMPLE','cage','genotype','sex','age'} <= set(data), 'Mouse identity fields are missing')
    require(not data.SAMPLE.duplicated().any() and data.cage.notna().all(), 'Invalid animal/cage identity')
    genotype = data.genotype.map({'WT': 0, '5X': 1})
    male = data.sex.str.upper().map({'F': 0, 'M': 1})
    require(genotype.notna().all() and male.notna().all(), 'Unknown genotype or sex')
    return pd.DataFrame({'sample': data.SAMPLE, 'cage': pd.factorize(data.cage, sort=False)[0],
                         'genotype': genotype.astype(int), 'age': pd.to_numeric(data.age),
                         'male': male.astype(int)})


def layouts(include_actual=True):
    answer = {}
    if include_actual:
        require(ACTUAL_INPUT is not None, 'Cage benchmark requires newly prepared --data-work')
        answer["actual_linear"] = observed_layout(ACTUAL_INPUT)
    for index, (size, mixed, fraction) in enumerate(itertools.product((2, 4, 6), (0., .5, 1.), (1/3, .5))):
        if size == 6 and fraction == .5:
            answer[f"n72_s6_m{int(mixed*100):03d}_g36"] = synthetic(size, mixed, fraction, 3252026+index)
    return answer


def cage_geometry(d, name):
    g = d.genotype.to_numpy(float)
    nuisance = [np.ones(len(d))]
    if "age" in d:
        nuisance += [d.age.to_numpy(float)-d.age.mean(), d.male.to_numpy(float)]
    z0 = np.column_stack(nuisance)
    x0 = np.column_stack([g, z0])
    zf = np.column_stack([z0, pd.get_dummies(d.cage, dtype=float).to_numpy()])
    means_matrix = np.column_stack([pd.Series(v).groupby(d.cage.to_numpy()).transform("mean") for v in [g, *nuisance[1:]]])
    zm = np.column_stack([z0, means_matrix])
    r0 = g-basis(z0)@(basis(z0).T@g)
    rf = g-basis(zf)@(basis(zf).T@g)
    rm = g-basis(zm)@(basis(zm).T@g)
    ss0, ssf, ssm = float(r0@r0), float(rf@rf), float(rm@rm)
    fixed_ok = ssf > 1e-9
    require(np.linalg.matrix_rank(x0) == x0.shape[1], "Baseline layout rank deficiency")
    require(abs(ssf-ssm) < 1e-8, "Mundlak and fixed-effect within information differs")
    if name == "actual_linear":
        require(len(d) == 74 and d.cage.nunique() == 36 and g.sum() == 38 and
                d.groupby("cage").genotype.nunique().eq(2).sum() == 17, "Observed housing manifest mismatch")
        require(abs(ss0-18.437395659432383) < 1e-8 and abs(ssf-8.305243445692884) < 1e-8,
                "Observed geometry differs from the manuscript")
    return dict(layout=name, n=len(d), cages=int(d.cage.nunique()),
                mixed_cages=int(d.groupby("cage").genotype.nunique().eq(2).sum()),
                SS_baseline=ss0, SS_fixed=ssf, SS_mundlak=ssm, fixed_estimable=fixed_ok,
                mundlak_estimable=bool(ssm > 1e-9),
                baseline_rank=int(np.linalg.matrix_rank(x0)),
                fixed_rank=int(np.linalg.matrix_rank(np.column_stack([g, zf]))),
                mundlak_rank=int(np.linalg.matrix_rank(np.column_stack([g, zm]))))


R_CODE = (Path(__file__).resolve().parents[1] / "analysis" / "benchmark.R").read_text(encoding="utf-8")


def r_setup(args, out):
    script = out / "runtime" / "analysis.R"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text(R_CODE, encoding="utf-8")
    exe = args.rscript or shutil.which("Rscript")
    require(exe and Path(exe).exists(), "Rscript not found; specify --rscript /path/to/Rscript")
    local = out / "runtime" / "R_library"
    local.mkdir(exist_ok=True)
    explicit = [Path(p).expanduser().resolve() for p in args.r_library]
    require(all(p.is_dir() for p in explicit), "An explicit --r-library is not a directory")
    found = []
    env = os.environ.copy()
    env["PAPER1_R_LIBRARIES"] = os.pathsep.join(map(str, [local, *explicit, *sorted(set(found))]))
    target, log = out / "runtime" / "R_environment.json", out / "runtime" / "R_check.log"
    print("Checking official R models; log:", log, flush=True)
    with log.open("w", encoding="utf-8") as stream:
        check = subprocess.run([str(exe), "--vanilla", str(script), "check", str(target),
                                "TRUE" if args.setup and not args.offline else "FALSE", str(local)],
                               env=env, stdout=stream, stderr=subprocess.STDOUT)
    if check.returncode:
        print(log.read_text(encoding="utf-8", errors="replace")[-8000:], flush=True)
        raise RuntimeError("R preflight failed; use --setup or an explicit --r-library. See " + str(log))
    # Exercise coefficient/variance/df interfaces before thousands of batches,
    # including the explicit-intercept Mundlak parameterization at balance.
    preflight = out / "runtime" / "model_preflight"
    preflight.mkdir(exist_ok=True)
    for label in ("actual_linear", "n72_s6_m100_g36"):
        units = preflight / (label + ".csv")
        layouts()[label].to_csv(units, index=False)
        result_path = preflight / (label + "_result.csv.gz")
        job = dict(id="preflight_"+label, spec=dict(scenario="preflight", layout=label,
                   rho_truth=0., shift=0., truth=0., seed=CAGE_SEED), start=0, end=1,
                   units=str(units), result=str(result_path), heartbeat=str(preflight / (label+".live.json")))
        job_path = preflight / (label + "_job.json")
        atomic_json(job_path, job)
        log_path = preflight / (label + ".log")
        with log_path.open("w", encoding="utf-8") as stream:
            trial = subprocess.run([str(exe), "--vanilla", str(script), "batch", str(job_path)],
                                   env=env, stdout=stream, stderr=subprocess.STDOUT)
        require(trial.returncode == 0, "R model preflight failed; inspect " + str(log_path))
        trial_rows = pd.read_csv(result_path)
        require(len(trial_rows) == 5 and trial_rows.status.eq("available").all(),
                "R model preflight returned unusable fits; inspect " + str(result_path))
    return dict(executable=str(Path(exe).resolve()), script=str(script), env=env, versions=read_json(target))


def receipt_ok(out, job, protocol_hash):
    path = Path(out) / "batches" / (job["id"] + ".csv.gz")
    receipt = read_json(str(path) + ".json")
    if not receipt:
        return False
    expected_rows = 12 if job["phase"] == "paired" else (job["end"]-job["start"])*len(METHODS)
    require(receipt.get("protocol_hash") == protocol_hash and receipt.get("job_id") == job["id"],
            "Receipt belongs to another protocol: " + str(path))
    require(path.is_file() and receipt.get("bytes") == path.stat().st_size and
            receipt.get("rows") == expected_rows and receipt.get("sha256") == sha(path),
            "Damaged batch; preserve it for inspection and use a new output or remove this batch AND its receipt: " + str(path))
    return True


def run_job(job, out, protocol_hash, r_info=None):
    start_time = time.monotonic()
    load_science()
    out = Path(out)
    path = out / "batches" / (job["id"] + ".csv.gz")
    heartbeat = out / "logs" / (job["id"] + ".live.json")
    atomic_json(heartbeat, dict(job_id=job["id"], phase=job["phase"], done=0,
                                total=job["end"]-job["start"], updated=time.time(), pid=os.getpid()))
    if job["phase"] == "paired":
        frame = paired_batch(job)
        atomic_csv(frame, path)
    else:
        config = dict(job, units=str(out / "inputs" / (job["spec"]["layout"] + ".csv")),
                      result=str(path), heartbeat=str(heartbeat))
        jp = out / "logs" / (job["id"] + ".job.json")
        atomic_json(jp, config)
        with (out / "logs" / (job["id"] + ".log")).open("w", encoding="utf-8") as stream:
            result = subprocess.run([r_info["executable"], "--vanilla", r_info["script"], "batch", str(jp)],
                                    env=r_info["env"], stdout=stream, stderr=subprocess.STDOUT)
        require(result.returncode == 0, "R batch failed; inspect logs/" + job["id"] + ".log")
        frame = pd.read_csv(path)
        require(len(frame) == (job["end"]-job["start"])*5 and
                not frame.duplicated(["replicate", "model"]).any(), "Incomplete R batch")
        require(set(frame.model) == set(METHODS) and
                frame.status.isin(["available", "failed", "not_estimable"]).all(), "Unknown R model/status")
    receipt = dict(job_id=job["id"], protocol_hash=protocol_hash, bytes=path.stat().st_size,
                   sha256=sha(path), rows=len(frame), repetitions=job["end"]-job["start"],
                   seconds=time.monotonic()-start_time, completed=time.time(), phase=job["phase"])
    atomic_json(str(path) + ".json", receipt)
    heartbeat.unlink(missing_ok=True)
    return receipt


def make_jobs(args):
    jobs, designs, cages = [], [], []
    if args.module in ("all", "paired"):
        grid = itertools.product((40, 80, 160), (0., .3, .7), ("ordinary", "high"), ("uniform", "fixed_gradient"))
        if args.smoke:
            grid = [(40, 0., "ordinary", "uniform"), (40, .3, "ordinary", "uniform"),
                    (40, .7, "high", "fixed_gradient")]
        for i, (n, r, leverage, weights) in enumerate(grid):
            base = dict(design_index=i, design=f"d{i:03d}", n=n, partial_correlation=r,
                        leverage=leverage, weights=weights)
            g = generic_design(base)
            designs.append(dict(**base, df=g["df"], max_leverage=float(g["h"].max()),
                                unit_variance=g["unit"], condition=float(np.linalg.cond(g["x1"]*np.sqrt(g["w"][:, None])))))
            if g["unit"] >= 1e-14:
                c = g["d"] @ g["x1"][:, -3:]
                index, v, tolerance = cancellation_direction(c)
                designs[-1].update(cancellation_coordinate=index, cancellation_tolerance=tolerance,
                                   cancellation_projection=json.dumps(c.tolist()),
                                   cancellation_direction=json.dumps(v.tolist()))
            for error in ("matched_gaussian", "heteroskedastic_gaussian", "matched_student_t5"):
                spec = dict(base, error_law=error)
                for lo in range(0, args.paired_reps, args.batch_size):
                    jobs.append(dict(id=f"p_{i:03d}_{error}_{lo:06d}", phase="paired", spec=spec,
                                     start=lo, end=min(lo+args.batch_size, args.paired_reps)))
    if args.module in ("all", "cage"):
        grid = list(itertools.product(layouts(), (0., .2), (0., .5), (0., .3)))
        if args.smoke:
            grid = [grid[0], grid[1], grid[6], grid[8], grid[16], grid[24], grid[30]]
        full_grid = list(itertools.product(layouts(), (0., .2), (0., .5), (0., .3)))
        for name, rho, shift, truth in grid:
            index = full_grid.index((name, rho, shift, truth))
            spec = dict(scenario=f"s{index:03d}", layout=name, rho_truth=rho, shift=shift,
                        truth=truth, seed=CAGE_SEED+index*100000)
            cages.append(spec)
            for lo in range(0, args.cage_reps, args.batch_size):
                jobs.append(dict(id=f"c_{index:03d}_{lo:06d}", phase="cage", spec=spec,
                                 start=lo, end=min(lo+args.batch_size, args.cage_reps)))
    return jobs, designs, cages


def from_sufficient(s):
    n, total = int(s["available"]), int(s["planned"])
    out = dict(planned=total, available=n, not_estimable=int(s["not_estimable"]),
               numerical_failed=int(s["numerical_failed"]), availability=n/total,
               coverage_all_planned=s["covered"]/total, rejection_all_planned=s["rejected"]/total)
    if n:
        bias = s["sum_error"]/n
        sd = math.sqrt(max(0., (s["sum_error2"]-n*bias*bias)/(n-1))) if n > 1 else float("nan")
        out.update(bias=bias, bias_MCSE=sd/math.sqrt(n), empirical_SD=sd, mean_SE=s["sum_SE"]/n,
                   RMS_SE=math.sqrt(s["sum_SE2"]/n), coverage=s["covered"]/n,
                   rejection=s["rejected"]/n, mean_width=s["sum_width"]/n)
        out["coverage_low"], out["coverage_high"] = wilson(s["covered"], n)
        out["rejection_low"], out["rejection_high"] = wilson(s["rejected"], n)
        out["rejection_MCSE"] = math.sqrt(out["rejection"]*(1-out["rejection"])/n)
    else:
        out.update({k: float("nan") for k in ("bias", "bias_MCSE", "empirical_SD", "mean_SE", "RMS_SE",
                   "coverage", "rejection", "mean_width", "coverage_low", "coverage_high",
                   "rejection_low", "rejection_high", "rejection_MCSE")})
    return out


def summarise(out, jobs, protocol):
    final = out / "final_results"
    final.mkdir(exist_ok=True)
    results = {}
    paired = [pd.read_csv(out / "batches" / (j["id"] + ".csv.gz")) for j in jobs if j["phase"] == "paired"]
    if paired:
        all_p = pd.concat(paired, ignore_index=True)
        keys = ["design", "n", "partial_correlation", "leverage", "weights", "error_law", "mean_case", "method", "status"]
        sums = ["planned", "available", "not_estimable", "numerical_failed", "sum_error", "sum_error2",
                "sum_SE", "sum_SE2", "covered", "rejected", "sum_width", "joint_F_rejected"]
        rows = []
        for values, part in all_p.groupby(keys, dropna=False, sort=True):
            s = part[sums].sum()
            row = dict(zip(keys, values))
            row.update(from_sufficient(s))
            row.update({k: part[k].iloc[0] for k in ("df", "max_leverage", "unit_variance", "true_delta", "oracle_SD")})
            row["joint_F_rejection"] = s.joint_F_rejected/s.planned if row["status"] == "regular" else float("nan")
            rows.append(row)
        tab = pd.DataFrame(rows)
        atomic_csv(tab, final / "paired_summary.csv")
        results["paired"] = tab
    cage_parts = [pd.read_csv(out / "batches" / (j["id"] + ".csv.gz")) for j in jobs if j["phase"] == "cage"]
    if cage_parts:
        c = pd.concat(cage_parts, ignore_index=True)
        require(not c.duplicated(["scenario", "replicate", "model"]).any(), "Repeated cage records")
        identity = c.mundlak_fixed_beta_difference.dropna()
        require(identity.empty or identity.abs().max() < 1e-6, "Mundlak coefficient identity failed")
        rows = []
        keys = ["scenario", "layout", "rho_truth", "shift", "truth", "model"]
        for values, part in c.groupby(keys, sort=True):
            valid = part.status.eq("available")
            d = part.loc[valid]
            covered = (d.low <= d.truth) & (d.high >= d.truth)
            s = dict(planned=len(part), available=len(d), not_estimable=int(part.status.eq("not_estimable").sum()),
                     numerical_failed=int(part.status.eq("failed").sum()), sum_error=float((d.beta-d.truth).sum()),
                     sum_error2=float(((d.beta-d.truth)**2).sum()), sum_SE=float(d.SE.sum()),
                     sum_SE2=float((d.SE**2).sum()), covered=int(covered.sum()), rejected=int(d.p.lt(.05).sum()),
                     sum_width=float((d.high-d.low).sum()))
            row = dict(zip(keys, values))
            row.update(from_sufficient(s))
            boundary = d.boundary.astype("string").str.lower().map({"true": 1., "false": 0.})
            row.update(rejection_interpretation="type_I_error" if row["truth"] == 0 else "power",
                       boundary_evaluable=int(boundary.notna().sum()),
                       boundary_rate=float(boundary.mean()) if boundary.notna().any() else float("nan"))
            rows.append(row)
        tab = pd.DataFrame(rows)
        for method in METHODS:
            selected = tab.loc[tab.model.eq(method)]
            require(selected.available.sum() > 0,
                    f"No usable {method} fit anywhere; inspect R logs before interpreting this run.")
        atomic_csv(tab, final / "cage_summary.csv")
        failures = c.loc[~c.status.eq("available") | c.warning.fillna("").ne("")]
        atomic_csv(failures, final / "cage_failures_and_warnings.csv.gz")
        atomic_csv(pd.DataFrame([dict(paired_available=len(identity), max_abs_beta_difference=float(identity.abs().max()) if len(identity) else float("nan"))]),
                   final / "mundlak_fixed_identity.csv")
        # Paired Monte Carlo differences use only common available repetitions;
        # missingness differences are reported separately, not silently matched away.
        pair_rows = []
        for scenario, part in c.groupby("scenario"):
            for comparator in METHODS[:-1]:
                a = part.loc[part.model.eq("mundlak_dream")].set_index("replicate")
                b = part.loc[part.model.eq(comparator)].set_index("replicate")
                require(a.index.equals(b.index), "Unmatched Monte Carlo repetitions")
                good = a.status.eq("available") & b.status.eq("available")
                row = dict(scenario=scenario, reference=comparator, planned=len(a), common_available=int(good.sum()),
                           availability_difference=float(a.status.eq("available").mean()-b.status.eq("available").mean()))
                for metric, va, vb in [("rejection", a.p.lt(.05), b.p.lt(.05)),
                                       ("coverage", a.low.le(a.truth)&a.high.ge(a.truth), b.low.le(b.truth)&b.high.ge(b.truth)),
                                       ("squared_error", (a.beta-a.truth)**2, (b.beta-b.truth)**2),
                                       ("width", a.high-a.low, b.high-b.low)]:
                    diff = va[good].astype(float)-vb[good].astype(float)
                    row[metric+"_difference"] = float(diff.mean()) if len(diff) else float("nan")
                    row[metric+"_difference_MCSE"] = float(diff.std(ddof=1)/np.sqrt(len(diff))) if len(diff)>1 else float("nan")
                pair_rows.append(row)
        atomic_csv(pd.DataFrame(pair_rows), final / "cage_paired_method_differences.csv")
        results["cage"] = tab
    for name in ("paired_geometry.csv", "cage_geometry.csv", "cage_scenarios.csv"):
        if (out / "inputs" / name).exists():
            shutil.copy2(out / "inputs" / name, final / name)
    shutil.copytree(out / "inputs", final / "inputs", dirs_exist_ok=True)
    atomic_json(final / "protocol.json", protocol)
    make_figures(final, results)
    numerical = sum(int(d.numerical_failed.sum()) for d in results.values())
    manifest = [{"path": str(p.relative_to(final)), "bytes": p.stat().st_size, "sha256": sha(p)}
                for p in sorted(final.rglob("*")) if p.is_file() and p.name != "file_manifest.csv"]
    atomic_csv(pd.DataFrame(manifest), final / "file_manifest.csv")
    return final, numerical


def make_figures(final, results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    colors = ["#555555", "#0072B2", "#009E73", "#CC79A7", "#D55E00"]
    if "paired" in results:
        d = results["paired"]
        d = d.loc[d.status.eq("regular") & d.mean_case.eq("quality_zero")]
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), constrained_layout=True)
        for ax, law in zip(axes, ("matched_gaussian", "heteroskedastic_gaussian", "matched_student_t5")):
            for i, method in enumerate(("model_based", "HC3", "independence_comparator")):
                z = d.loc[d.error_law.eq(law) & d.method.eq(method)]
                ax.scatter(z.max_leverage, z.rejection, s=16, alpha=.7, label=method, color=colors[i])
            ax.axhline(.05, color="black", ls="--", lw=.8)
            ax.set(xlabel="Maximum leverage", ylabel="Pointwise null rejection", title=law.replace("_", " "))
        axes[0].legend(fontsize=7)
        for ext in ("pdf", "png"):
            fig.savefig(final / ("paired_calibration."+ext), dpi=300)
        plt.close(fig)
    if "cage" in results:
        d = results["cage"]
        d = d.loc[d.layout.eq("actual_linear") & d.truth.eq(0)]
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), constrained_layout=True)
        cases = d[["rho_truth", "shift"]].drop_duplicates().sort_values(["rho_truth", "shift"])
        labels = [f"rho={r:g}, shift={s:g}" for r, s in cases.itertuples(index=False, name=None)]
        for i, method in enumerate(METHODS):
            z = d.loc[d.model.eq(method)].sort_values(["rho_truth", "shift"])
            xx = np.arange(len(z))+(i-2)*.13
            for ax, metric in zip(axes, ("bias", "rejection", "mean_width")):
                if metric == "rejection":
                    err = np.vstack([z.rejection-z.rejection_low, z.rejection_high-z.rejection])
                elif metric == "bias":
                    err = 1.96*z.bias_MCSE
                else:
                    err = None
                ax.errorbar(xx, z[metric], yerr=err, fmt="o", ms=3, capsize=2, color=colors[i], label=method)
                ax.set_xticks(np.arange(len(labels)), labels, rotation=35, ha="right")
                ax.set_ylabel(metric.replace("_", " "))
        axes[0].axhline(0, ls="--", lw=.8, color="black")
        axes[1].axhline(.05, ls="--", lw=.8, color="black")
        axes[0].legend(fontsize=7)
        for ext in ("pdf", "png"):
            fig.savefig(final / ("cage_model_comparison."+ext), dpi=300)
        plt.close(fig)


def status(out):
    protocol = read_json(out / "protocol.json")
    if not protocol:
        print("No protocol yet. Check the launch terminal and runtime/R_check.log.")
        return
    print("Recorded state:", read_json(out / "run_state.json").get("status", "unknown"))
    with output_lock(out, probe=True) as free:
        print("Writer lock:", "INACTIVE" if free else "ACTIVE")
    for phase in ("paired", "cage"):
        jobs = [j for j in protocol["jobs"] if j["phase"] == phase]
        if not jobs:
            continue
        done, damaged, seconds = 0, 0, []
        for j in jobs:
            try:
                if receipt_ok(out, j, protocol["hash"]):
                    done += 1
                    seconds.append(read_json(out / "batches" / (j["id"]+".csv.gz.json"))["seconds"])
            except ValueError:
                damaged += 1
        print(f"{phase}: {done}/{len(jobs)} receipt/hash-checked batches ({100*done/len(jobs):.1f}%), damaged={damaged}")
        state = read_json(out / "run_state.json")
        workers = state.get("r_workers" if phase == "cage" else "workers", 1)
        if state.get("status") == "RUNNING" and len(seconds) >= min(workers, 5):
            median = sorted(seconds)[len(seconds)//2]
            print(f"  Rough phase-only remaining: {(len(jobs)-done)*median/workers/3600:.2f} h at {workers} workers; excludes other phases/setup/summaries.")
    live = [read_json(p) for p in (out / "logs").glob("*.live.json")]
    live = [d for d in live if time.time()-d.get("updated", 0) < 120]
    print("Recent batch heartbeats:", len(live), "(not a process count)")
    print("Compact results:", out / "final_results")


def self_test():
    load_science()
    g = generic_design(dict(n=80, partial_correlation=.4, leverage="ordinary", weights="fixed_gradient", design_index=123))
    rng = np.random.default_rng(31)
    y = rng.normal(size=(80, 5))
    b0 = np.linalg.lstsq(g["x0"]*np.sqrt(g["w"][:, None]), y*np.sqrt(g["w"][:, None]), rcond=None)[0]
    b1 = np.linalg.lstsq(g["x1"]*np.sqrt(g["w"][:, None]), y*np.sqrt(g["w"][:, None]), rcond=None)[0]
    np.testing.assert_allclose(g["d"]@y, b1[1]-b0[1], atol=1e-12)
    c = np.eye(g["x1"].shape[1])[1]-g["a0"][1]@g["x1"]
    np.testing.assert_allclose(c@np.linalg.inv(g["x1"].T@(g["w"][:, None]*g["x1"]))@c, g["unit"], atol=1e-12)
    for _, mu, truth in means(g):
        np.testing.assert_allclose(g["d"]@mu, truth, atol=1e-12)
    for name, d in layouts(include_actual=ACTUAL_INPUT is not None).items():
        z = cage_geometry(d, name)
        require(z["fixed_estimable"] == z["mundlak_estimable"], "Estimability test failed")
    spec = dict(design="test", design_index=18, n=40, partial_correlation=.3,
                leverage="high", weights="fixed_gradient", error_law="heteroskedastic_gaussian")
    one = paired_batch(dict(spec=spec, start=0, end=21))
    two = [paired_batch(dict(spec=spec, start=a, end=b)) for a, b in [(0, 8), (8, 21)]]
    for column in ("sum_error", "sum_error2", "sum_SE", "sum_SE2", "rejected", "covered"):
        np.testing.assert_allclose(one[column], two[0][column]+two[1][column], rtol=1e-11, atol=1e-11)
    print("PASS: WLS maps, paired covariance, constructed truths, synthetic Mundlak geometry, batch-invariant streams.")


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--project-root", type=Path, default=Path.cwd(), help="Only used for local R-library discovery")
    p.add_argument("--output", type=Path, default=None)
    p.add_argument("--data-work", type=Path, help="Fresh core output; observed allocation is read from mouse_units.csv")
    p.add_argument("--resume", action="store_true", help="Continue this same raw-input experiment only")
    p.add_argument("--workers", type=int, default=min(36, os.cpu_count() or 1), help="Python paired-grid processes")
    p.add_argument("--r-workers", type=int, default=None, help="Concurrent single-thread R batches; default min(workers,24)")
    p.add_argument("--paired-reps", type=int, default=5000)
    p.add_argument("--cage-reps", type=int, default=5000)
    p.add_argument("--batch-size", type=int, default=250)
    p.add_argument("--rscript", default=None)
    p.add_argument("--r-library", action="append", default=[], help="Repeat for explicit existing R libraries")
    p.add_argument("--module", choices=("all", "paired", "cage"), default="all")
    p.add_argument("--setup", action="store_true", help="Install missing dependencies only; no research-data downloads")
    p.add_argument("--offline", action="store_true", help="Forbid dependency installation")
    p.add_argument("--smoke", action="store_true", help="Small diagnostic run; NOT manuscript evidence")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--status", action="store_true", help="Read-only progress/ETA check")
    return p


def main():
    global ACTUAL_INPUT
    args = parser().parse_args()
    args.project_root = args.project_root.expanduser().resolve()
    args.output = (args.output or args.project_root.parent / "paper1_benchmark_results").expanduser().resolve()
    if args.status:
        status(args.output)
        return 0
    require(args.workers > 0 and args.batch_size > 0 and args.paired_reps > 1 and args.cage_reps > 1, "Invalid workers/repetition/batch setting")
    require(args.cage_reps <= 100000, "Cage scenario seeds are spaced by 100000; larger repetition counts would overlap streams")
    args.r_workers = args.r_workers if args.r_workers is not None else min(args.workers, 24)
    require(args.r_workers > 0, "Invalid R worker count")
    bootstrap(args)
    load_science()
    if args.self_test:
        self_test()
        return 0
    if args.smoke:
        args.paired_reps, args.cage_reps, args.batch_size = 30, 4, 15
    from pipeline.lifecycle import claim, current_source, parent_run, code_identity
    source_id = parent_run()
    if args.module != "paired":
        require(args.data_work is not None, 'Specify the freshly generated --data-work for the cage benchmark')
        require(source_id is not None and args.data_work.resolve() == args.output.parent / 'core',
                'Use run.py: cage inputs must be prepared in this same composite experiment.')
        source_id = current_source(args.data_work, source_id)
        ACTUAL_INPUT = args.data_work.resolve() / 'prepared' / 'mouse_units.csv'
        require(ACTUAL_INPUT.is_file(), 'Fresh mouse identity table is missing')
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    with output_lock(out):
        run_id = claim(out, 'benchmark', args.resume, source_id,
                       code_identity(Path(__file__).resolve().parents[1]))
        for name in ("inputs", "batches", "logs", "runtime"):
            (out / name).mkdir(exist_ok=True)
        r_info = r_setup(args, out) if args.module != "paired" else None
        import scipy
        import matplotlib
        scientific = dict(version=VERSION, run_id=run_id, source_sha256=sha(Path(__file__)),
                          r_source_sha256=sha(Path(__file__).resolve().parents[1] / 'analysis' / 'benchmark.R'), module=args.module,
                          paired_reps=args.paired_reps, cage_reps=args.cage_reps, batch_size=args.batch_size,
                          paired_seed=PAIRED_SEED, cage_seed=CAGE_SEED, smoke=args.smoke,
                          cancellation_direction_rule=CANCELLATION_RULE,
                          metadata_sha256=sha(ACTUAL_INPUT) if ACTUAL_INPUT else None,
                          metadata_source=str(ACTUAL_INPUT) if ACTUAL_INPUT else 'synthetic paired design only',
                          paired_rng="SeedSequence [seed,1000,design] for designs; [seed,2000,design,error,repetition] for errors",
                          python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                          pandas=pd.__version__, matplotlib=matplotlib.__version__,
                          R=None if not r_info else r_info["versions"],
                          analysis_scale="fixed-design normalized-response; no count, voom, eBayes or BH recalibration",
                          mundlak="genotype plus all baseline covariate cage means; target must remain estimable",
                          boundary="dream tau/(tau+sigma2)<1e-8; valid boundary fits retained")
        protocol_hash = digest(scientific)
        existing = read_json(out / "protocol.json")
        require(not existing or existing.get("hash") == protocol_hash,
                "Protocol/code/environment changed. Use a NEW --output; worker counts alone may change when resuming.")
        jobs, designs, cages = make_jobs(args)
        active_layouts = layouts() if args.module != 'paired' else {}
        for name, units in active_layouts.items():
            units.to_csv(out / "inputs" / (name + ".csv"), index=False)
        if active_layouts:
            pd.DataFrame([cage_geometry(d, name) for name, d in active_layouts.items()]).to_csv(out / "inputs" / "cage_geometry.csv", index=False)
        if designs:
            pd.DataFrame(designs).to_csv(out / "inputs" / "paired_geometry.csv", index=False)
        if cages:
            pd.DataFrame(cages).to_csv(out / "inputs" / "cage_scenarios.csv", index=False)
        protocol = dict(hash=protocol_hash, scientific=scientific, jobs=jobs)
        if not existing:
            atomic_json(out / "protocol.json", protocol)
        started = time.time()
        state = dict(status="RUNNING", pid=os.getpid(), workers=args.workers, r_workers=args.r_workers,
                     started=started, protocol_hash=protocol_hash)
        atomic_json(out / "run_state.json", state)
        try:
            for phase in ("paired", "cage"):
                tasks = [j for j in jobs if j["phase"] == phase]
                pending = [j for j in tasks if not receipt_ok(out, j, protocol_hash)]
                workers = args.workers if phase == "paired" else args.r_workers
                print(f"{phase}: verified {len(tasks)-len(pending)}/{len(tasks)}; pending {len(pending)}, workers={workers}", flush=True)
                state.update(phase=phase, updated=time.time())
                atomic_json(out / "run_state.json", state)
                if not pending:
                    continue
                if workers == 1:
                    for index, job in enumerate(pending, len(tasks)-len(pending)+1):
                        receipt = run_job(job, str(out), protocol_hash, r_info)
                        print(f"{phase}: {index}/{len(tasks)} batch={job['id']} seconds={receipt['seconds']:.1f}", flush=True)
                    continue
                executor_type = cf.ProcessPoolExecutor if phase == "paired" else cf.ThreadPoolExecutor
                # Bound dispatch so interruption does not leave thousands of queued jobs.
                pool = executor_type(max_workers=workers)
                running = {}
                iterator = iter(pending)
                completed = len(tasks)-len(pending)
                try:
                    for _ in range(min(len(pending), workers*2)):
                        j = next(iterator)
                        running[pool.submit(run_job, j, str(out), protocol_hash, r_info)] = j
                    while running:
                        done, _ = cf.wait(running, timeout=30, return_when=cf.FIRST_COMPLETED)
                        if not done:
                            print(f"{phase}: completed {completed}/{len(tasks)}; waiting for active batches", flush=True)
                        for future in done:
                            j = running.pop(future)
                            receipt = future.result()
                            completed += 1
                            print(f"{phase}: {completed}/{len(tasks)} batch={j['id']} seconds={receipt['seconds']:.1f}", flush=True)
                            nxt = next(iterator, None)
                            if nxt is not None:
                                running[pool.submit(run_job, nxt, str(out), protocol_hash, r_info)] = nxt
                finally:
                    # Let only in-flight batches finish and publish complete receipts.
                    for f in running:
                        f.cancel()
                    pool.shutdown(wait=True, cancel_futures=True)
            require(all(receipt_ok(out, j, protocol_hash) for j in jobs), "Some planned batches are missing")
            require(sha(Path(__file__)) == scientific["source_sha256"] and
                    sha(Path(__file__).resolve().parents[1] / "analysis" / "benchmark.R") == scientific["r_source_sha256"],
                    "Source code changed during this run. Preserve batches and restore the original code before resuming.")
            final, failures = summarise(out, jobs, protocol)
            state.update(status="COMPLETED_WITH_FIT_FAILURES" if failures else "COMPLETED",
                         updated=time.time(), seconds=time.time()-started, final_results=str(final))
            atomic_json(out / "run_state.json", state)
            print("FINAL RESULTS", final, flush=True)
            return 0
        except BaseException as exc:
            state.update(status="INTERRUPTED" if isinstance(exc, KeyboardInterrupt) else "FAILED",
                         updated=time.time(), error=str(exc))
            atomic_json(out / "run_state.json", state)
            raise


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted. Add --resume with identical settings to continue this same run.", file=sys.stderr)
        raise SystemExit(130)
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
