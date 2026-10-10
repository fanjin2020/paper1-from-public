# Targeted Review of dream Convergence Warnings

[English run guide](../README_EN.md) | [Chinese version](DREAM_REPAIR.md)

## When to Use This Entry Point

Version 1.2.2 replaces only raw-fit rows that pass a two-optimizer agreement check. Unaffected genes retain their original raw estimates. Negative-eigenvalue and convergence warnings are localized and evaluated rather than suppressed.

The targeted entry point reads affected tasks from the original diagnostics; it does not hard-code a gene to repair. It currently accepts a **complete v1.2.0 work directory**, not a review ZIP without count matrices.

**For a new analysis, use `run.py`.** It already includes the same convergence-review logic. A completed v1.2.2 analysis does not need the v1.2.0 migration entry point.

## Review Procedure

1. Verify SHA-256 records for original results, counts, and animal metadata. Use separate source and destination directories; retain the original files.
2. Verify R dependencies and versions, then run regression tests and the official dream API smoke test before using real counts.
3. Reconstruct affected inputs with the same background, sample order, design, TMM normalization, and voom settings.
4. Reproduce the default fit. Available estimates, missing positions, and statuses must agree with the original output within the specified tolerance.
5. Retry localized genes with both `bobyqa` and `nloptwrap`. Neither retry may retain unresolved warnings or background errors. Compare target coefficients, raw SEs, degrees of freedom, likelihood, residual scale, and cage variance.
6. For accepted genes, use the prespecified bobyqa result to replace only the raw-fit row, including its required coefficients, scale, degrees of freedom, covariance, and variance components. Validate object structure, design, and gene identity during the merge. Unknown fields, already moderated objects, and unmappable rows stop the merge.
7. Apply empirical Bayes moderation and interval extraction over the entire usable background, then apply BH correction over the original planned family. Do not insert a separately moderated single-gene p value into the old table.
8. On the same final raw fit, calculate an additional result excluding the restored gene to quantify its influence on the other genes' moderated SEs and p values.

A gene is not restored if optimizers disagree, warnings remain, or the raw objects cannot be merged safely. The procedure keeps the formula, cage term, convergence criteria, and optimizer-selection rule unchanged.

## Numerical Acceptance Rules

Optimizer agreement uses:

```text
abs(a - b) <= tolerance * (1 + max(abs(a), abs(b)))
```

| Quantity | Tolerance |
|---|---:|
| Coefficient | `1e-5` |
| Unmoderated SE | `1e-4` |
| Degrees of freedom | `1e-3` |
| Likelihood | `1e-6` |
| Residual scale | `1e-4` |
| Cage variance | `1e-4` |

These tolerances describe numerical agreement, not biological equivalence. Zero cage variance is allowed. Nonfinite estimates, nonpositive SEs, or nonpositive degrees of freedom are rejected.

Default-fit reproduction requires:

```text
max(abs(old - new) / (1 + abs(old))) < 1e-6
```

Missing positions and statuses are checked as well. The original R package versions must match so that a dependency upgrade is not mixed into the targeted comparison.

## Run the Targeted Review

Activate the verified environment and run from the new code directory containing `run_dream_repair.py`. Replace all example paths with real locations.

```bash
conda activate paper1
set -euo pipefail

python -u run_dream_repair.py \
  --source-work /path/to/v1.2.0/work_reanalysis_v1_2_0 \
  --work /path/to/new/work_dream_repair_v1_2_2 \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --library /path/to/existing/R_library \
  --workers 1 \
  2>&1 | tee dream_repair.log
```

This does not download data again. The destination uses symbolic links to access the two source mouse count matrices without rewriting them; keep the source directory available. Result tables are copied before being updated.

The targeted entry point does not support resumption after interruption. Use a new destination directory for another attempt. This differs from `run.py`, which supports `--resume` for the complete workflow.

## Outputs

| File, relative to the destination work directory | Purpose |
|---|---|
| `results/mouse_dream_retry_attempts.csv` | Raw optimizer estimates, warnings, and agreement checks |
| `results/dream_repair_default_reproduction.csv` | Default-fit reproduction against v1.2.0 |
| `results/dream_repair_before_after.csv.gz` | Before/after results for all background genes in affected models |
| `results/mouse_dream_background_impact.csv.gz` | Common moderation with and without restored genes; empty if none was restored |
| `results/mouse_dream_diagnostics.csv` | Final localization and recovery status |
| `results/mouse_fit_failures.csv` | Remaining unavailable fits |
| `repair_changes.json` | Updated and unchanged source-file records |
| `review_bundle.zip` | Results, logs, and executed code for review |

`accepted` means that the two optimizer solutions passed numerical checks. Final restoration requires `recovered_optimizer_agreement` in the diagnostics and `available` in the complete results.

The original `mouse_variance.csv.gz` comes from a separate variance-component fit and is not rewritten by this entry point. Cage variance from the retry itself is recorded in `mouse_dream_retry_attempts.csv`.

## Exit Codes

| Code | Meaning |
|---|---|
| `0` | Integrity checks passed and the reviewed fit positions are available |
| `3` | Computation ended with remaining unavailable fits or unresolved warnings; inspect diagnostics |
| `2` | Execution, environment, or source-reproduction checks failed |

Human models, no-cage/fixed/common-correlation mouse models, SE decomposition, simulations, and paired-test calibration are reused with checksum checks. The corrected exact-SE calculation is retained; the historical SE discrepancy and convergence review are separate issues.

## Implementation References

- [lme4 convergence guidance](https://lme4.github.io/lme4/reference/convergence.html): compare optimizers and stopping accuracy instead of suppressing warnings.
- [Official dream interface](https://gabrielhoffman.github.io/variancePartition/reference/dream-method.html): pass optimizer control through the public API while retaining REML and Satterthwaite degrees of freedom.

See [VALIDATION_EN.md](VALIDATION_EN.md) for release evidence. Each server run executes its own preflight; convergence of real-data fits is determined from that run's diagnostics.
