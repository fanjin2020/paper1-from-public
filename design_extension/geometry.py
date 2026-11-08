"""Outcome-free geometry; covariance-dependent results require supplied assumptions."""
import numpy as np


def basis(matrix):
    matrix = np.asarray(matrix, float)
    if matrix.ndim != 2 or not np.isfinite(matrix).all():
        raise ValueError("Nuisance matrix must be finite and two-dimensional")
    if not matrix.shape[1]:
        return np.zeros((len(matrix), 0))
    norm = np.linalg.norm(matrix, axis=0)
    matrix = matrix[:, norm > 0] / norm[norm > 0]
    if not matrix.shape[1]:
        return np.zeros((len(matrix), 0))
    u, s, _ = np.linalg.svd(matrix, full_matrices=False)
    return u[:, s > max(matrix.shape) * np.finfo(float).eps * s[0] * 100]


def residual(x, nuisance):
    q = basis(nuisance)
    return x - q @ (q.T @ x)


def components(x, groups, covariates=None):
    x = np.asarray(x, float)
    groups = np.asarray(groups)
    if x.ndim != 1 or len(x)<2 or groups.ndim!=1 or len(x) != len(groups) or not np.isfinite(x).all():
        raise ValueError("Invalid exposure/group alignment")
    if any(g is None or str(g).strip().lower() in ('','nan','none') for g in groups):
        raise ValueError("Missing group ID")
    _, ix = np.unique(groups, return_inverse=True)
    z = np.eye(ix.max() + 1)[ix]
    c = np.ones((len(x), 1))
    if covariates is not None:
        covariates = np.asarray(covariates, float)
        if covariates.ndim != 2 or len(covariates) != len(x):
            raise ValueError("Covariate dimensions differ")
        c = np.column_stack([c, covariates])
    return x, z, basis(c), ix


def design_diagnostics(x, groups, covariates=None):
    x, z, c, ix = components(x, groups, covariates)
    within = residual(x, z)
    base = residual(x, c)
    fixed = residual(x, np.column_stack([c, z]))
    total = float(np.sum((x - x.mean()) ** 2))
    ib, iw = float(base @ base), float(fixed @ fixed)
    threshold = 1e-10 * max(total, 1.)
    information_by_group = np.bincount(ix, weights=fixed ** 2)
    positive = information_by_group > threshold
    effective = (iw ** 2 / np.sum(information_by_group ** 2)) if iw > threshold else 0.
    df = len(x) - basis(np.column_stack([c, z, x])).shape[1]
    estimable = iw > threshold
    return dict(n=len(x), groups=z.shape[1], total_information=total,
                ewi=float(np.clip(within @ within / total,0.,1.)) if total > threshold else None,
                conditional_ewi=float(np.clip(iw / ib,0.,1.)) if ib > threshold else None,
                baseline_information=ib, within_information=iw,
                baseline_estimable=ib > threshold, within_estimable=estimable,
                fixed_residual_df=df, groups_with_conditional_information=int(positive.sum()),
                effective_information_groups=float(effective),
                largest_group_information_share=float(information_by_group.max() / iw) if estimable else None,
                homoskedastic_se_multiplier=float(np.sqrt(ib / iw)) if estimable else None,
                design_status=("target_not_estimable" if ib <= threshold else
                               "within_effect_not_estimable" if not estimable else
                               "insufficient_residual_df" if df < 5 else "within_effect_estimable"))


def linear_operator(x, groups, covariates=None, model="baseline", cage_sd=0.):
    """Known-covariance Gaussian sensitivity, NOT estimated dream inference."""
    if cage_sd < 0 or not np.isfinite(cage_sd):
        raise ValueError("cage_sd must be finite and nonnegative")
    x, z, c, _ = components(x, groups, covariates)
    means = z @ (z.T @ np.column_stack([c, x]) / z.sum(axis=0)[:, None])
    nuisance = c
    if model == "fixed":
        nuisance = np.column_stack([c, z])
    elif model == "mundlak_known_covariance":
        nuisance = np.column_stack([c, means])
    elif model not in ("baseline", "re_known_covariance"):
        raise ValueError("Unknown Gaussian model")
    if np.linalg.norm(residual(x, nuisance)) ** 2 < 1e-10 * max(1., np.sum((x-x.mean())**2)):
        return None
    covariance = np.eye(len(x))
    if model.endswith("known_covariance"):
        covariance = covariance + cage_sd ** 2 * (z @ z.T)
    inv = np.linalg.solve(covariance, np.eye(len(x)))
    design = np.column_stack([basis(nuisance), x])
    a = np.linalg.solve(design.T @ inv @ design, design.T @ inv)
    df = len(x) - design.shape[1]
    if df < 5:
        return None
    target = a[-1]
    g = means[:, -1] - x.mean()
    return dict(a=target, design=design, all_a=a, inv=inv, df=df,
                assumed_se_unit=float(np.sqrt(target @ covariance @ target)),
                gamma_bias_multiplier=float(target @ g),
                group_shift_l2_multiplier=float(np.linalg.norm(z.T @ target)))
