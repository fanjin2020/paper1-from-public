"""Prespecified design grid and independent, paired simulation streams."""
import csv
from pathlib import Path
import numpy as np
from scipy.stats import t
from .geometry import design_diagnostics, linear_operator

MODELS = ("baseline", "fixed", "re_known_covariance", "mundlak_known_covariance")


def read_design(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows or len({r["animal_id"] for r in rows}) != len(rows):
        raise ValueError("Empty metadata or duplicate animal ID")
    x = np.array([{"WT": 0., "5X": 1., "0": 0., "1": 1.}[r["genotype"]] for r in rows])
    groups = np.array([r["cage"] for r in rows])
    cov = np.array([[float(r["age"]), {"F": 0., "M": 1.}[r["sex"]]] for r in rows])
    return rows, x, groups, cov


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def create_designs(folder, seed):
    designs = []
    # Holding out ALL 18-cage designs prevents random-seed-only validation.
    for cages in (12, 18, 36):
        size = 72 // cages
        for level in range(4):
            mixed = cages * level // 3
            pure = cages - mixed
            if pure % 2 or size % 2:
                raise ValueError("Layout needs balanced even group sizes")
            counts = [size // 2] * mixed + [0] * (pure // 2) + [size] * (pure // 2)
            identifier = f"k{cages:02d}_w{level}"
            rng = np.random.default_rng(np.random.SeedSequence([seed, cages, level, 41]))
            rows = []
            for k, mutants in enumerate(counts):
                genotype = np.array([1] * mutants + [0] * (size - mutants))
                rng.shuffle(genotype)
                for v in genotype:
                    i = len(rows)
                    rows.append(dict(animal_id=f"a{i:03d}", cage=f"c{k:03d}",
                                     genotype="5X" if v else "WT", age=float(rng.uniform(3, 12)),
                                     sex="M" if rng.integers(2) else "F", slide=f"s{i//12}"))
            path = folder / (identifier + ".csv")
            write_csv(path, rows)
            _, x, g, cov = read_design(path)
            designs.append(dict(design=identifier, split="holdout" if cages == 18 else "training",
                                path=str(path), **design_diagnostics(x, g, cov)))
    return designs


def scenarios(profile):
    rows = []
    for pi in ([.1] if profile == "pilot" else [0., .1]):
        for sd in (0., .3):
            for gamma in (0., .5):
                rows.append(dict(scenario=len(rows), cage_log2_sd=sd, cage_mean_amplitude=gamma,
                                 de_fraction=pi, effect_log2=.5, mean_affected_fraction=.2))
    return rows


def bh(p):
    p = np.asarray(p, float)
    available = np.isfinite(p)
    q = np.ones(len(p))
    order = np.argsort(np.where(available, p, 1.))
    values = np.minimum.accumulate((np.where(available, p, 1.)[order] * len(p) /
                                   np.arange(1, len(p)+1))[::-1])[::-1]
    q[order] = np.minimum(1., values)
    q[~available] = np.nan
    return q


def gaussian_job(job, output):
    _, x, groups, cov = read_design(job["design_path"])
    from .geometry import components
    _, z, _, _ = components(x, groups, cov)
    n = len(x)
    s = job["scenario"]
    ops = {m: linear_operator(x, groups, cov, m, s["cage_log2_sd"]) for m in MODELS}
    g = z @ (z.T @ x / z.sum(axis=0)) - x.mean()
    rows, sensitivity = [], []
    for m, op in ops.items():
        sensitivity.append(dict(design=job["design"], scenario=s["scenario"], model=m,
            estimable=op is not None,
            gamma_bias_multiplier=op["gamma_bias_multiplier"] if op else None,
            group_shift_l2_multiplier=op["group_shift_l2_multiplier"] if op else None,
            assumed_se_unit=op["assumed_se_unit"] if op else None))
    for rep in range(job["start"], job["end"]):
        # Same underlying normals/support across gamma and covariance conditions.
        rng = np.random.default_rng(np.random.SeedSequence([job["seed"], job["design_index"], rep, 921]))
        ng = job["genes"]
        order = rng.permutation(ng)
        beta = np.zeros(ng)
        alt = order[:int(round(ng*s["de_fraction"]))]
        beta[alt] = np.resize(np.array([-1., 1.]), len(alt))*s['effect_log2']
        affected = rng.permutation(ng)[:max(1, int(s['mean_affected_fraction']*ng))]
        gamma = np.zeros(ng)
        gamma[affected] = np.resize(np.array([-1., 1.]), len(affected))*s["cage_mean_amplitude"]
        u = rng.normal(size=(z.shape[1], ng)) * s["cage_log2_sd"]
        y = np.outer(x, beta) + np.outer(g, gamma) + z @ u + rng.normal(size=(n, ng))
        for model, op in ops.items():
            row = dict(design=job["design"], split=job["split"], scenario=s["scenario"],
                       replicate=rep, model=model, available=0., FDP=None, power=None,
                       coverage=None, mean_bias=None, bias_prediction_rmse=None)
            if op:
                estimated = op["a"] @ y
                e = y - op["design"] @ (op["all_a"] @ y)
                scale = np.sum(e*(op["inv"]@e), axis=0)/op["df"]
                se = np.sqrt(np.maximum(scale, 0.))*op["assumed_se_unit"]
                p = 2*t.sf(np.abs(estimated/se), op["df"])
                discovered = bh(p) < .05
                alternative = beta != 0
                error = estimated-beta
                predicted = gamma*op["gamma_bias_multiplier"]
                half = t.ppf(.975, op["df"])*se
                row.update(available=1., FDP=float(np.sum(discovered & ~alternative)/max(1, discovered.sum())),
                           power=float(discovered[alternative].mean()) if alternative.any() else None,
                           coverage=float((np.abs(error)<=half).mean()), mean_bias=float(error.mean()),
                           bias_prediction_rmse=float(np.sqrt(np.mean((error-predicted)**2))),
                           bias_at_positive_gamma=float(error[gamma>0].mean()) if (gamma>0).any() else None,
                           predicted_bias_at_positive_gamma=float(predicted[gamma>0].mean()) if (gamma>0).any() else None)
            rows.append(row)
    write_csv(output / "metrics.csv", rows)
    write_csv(output / "sensitivity.csv", sensitivity)
