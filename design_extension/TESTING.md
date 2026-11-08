# Validation Record

Date: 2026-11-01. Extension version: 0.1.0.

## Completed Locally

- 17 Python unit tests passed using Python 3.12.2, NumPy 2.2.6 and SciPy 1.15.3.
- Tested within-information identities, pure/balanced boundaries, covariate
  collinearity, coefficient bias formulas, nuisance group-shift bounds,
  Mundlak/FE coefficient equivalence under the stated known covariance,
  row/group relabeling and covariate rescaling.
- Tested held-out-result exclusion, pilot abstention, paired mean-shift
  summaries, nonzero uncertainty after zero FDP events, planned-family BH,
  complete-arm failure detection, checksum corruption detection and protected
  output locations.
- Final standalone-package Gaussian smoke run completed 48 jobs / 48 datasets
  and generated its compact review ZIP. Resume with a different worker count
  found zero pending batches and performed no new fits.
- R 4.4.2 passed source parsing, rank/estimability checks, deterministic negative
  binomial generation, truth support sizes, planned-family BH and unavailable
  inference accounting. These tests used a small constructed fixture, not
  research observations.
- All five frozen R engine files match their repository originals by SHA256.
- No software packages were installed for these checks. Shell entry uses LF
  line endings for Linux.

## Not Yet Validated on the Target Server

The local R installation lacks loadable edgeR, variancePartition, BiocParallel,
lme4 and lmerTest. The complete count-calibration and five-model fitting path
has therefore NOT been executed locally. Its source uses the existing
count-generation and fitting functions without modification, but the new grid
and orchestration still need the server pilot.

No formal count benchmark, held-out performance result, publication figure or
manuscript conclusion has been produced by this code delivery. The pilot is an
execution check; full-profile results require separate numerical review.

## Re-run Tests

From the package/repository root in an existing numerical environment:

```bash
python -m unittest discover -s design_extension/tests -v
Rscript --vanilla design_extension/tests/test_engine.R
```

The R fixture tests do not require dream. The actual count experiment does.
