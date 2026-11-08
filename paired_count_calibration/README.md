# End-to-end count calibration of paired coefficient differences

This extension refits the complete pipeline in every independent dataset:

`NB counts -> filterByExpr -> TMM -> voom -> base/full lmFit -> paired contrast -> BH`.

It calls the original `analysis/matched/functions.R::paired_weights` unchanged.
The paired model-based variance is unmoderated, as in that implementation;
the second test uses the full-model residuals and HC3 leverage correction.
Empirical-Bayes moderation is not applied to the paired difference variance.
Both shared estimated voom weights and design-specific estimated weights are tested.
The two coefficient estimates use the same logCPM response and retained genes.
Each dataset checks that their difference agrees with the original paired kernel.

## Run on the existing server environment

No downloads, package installation, old fitted-result imports or Git operations.
Python uses only its standard library. R needs edgeR, limma and jsonlite.
The default version check uses the existing manuscript R/edgeR/limma pins.
An explicit `--allow-version-drift` records a different environment instead of silently accepting it.
Do not mix environments within a run.

```bash
conda activate paper1
bash run_paired_count_calibration.sh --profile full --workers 36 \
  --output "$HOME/adjscope/p1_paired_count_full"
```

Full: 24 scenarios x 1,000 independent datasets x 3,000 genes.
Pilot: 24 x 20 x 1,000. Smoke: four scenarios x one dataset x 250 genes.
These are independent simulations; the full run does not reuse pilot outcomes.
`--reps` and `--genes` override these counts and are recorded in the protocol.
For a short environment test, use `--profile smoke` and a different output directory.

```bash
# Check cumulative receipt-verified DATASETS, not batches since the last restart.
bash run_paired_count_calibration.sh --status \
  --output "$HOME/adjscope/p1_paired_count_full"

# Resume after Ctrl-C or a server restart. Use the SAME actual output directory.
bash run_paired_count_calibration.sh --resume --workers 36 \
  --output "$HOME/adjscope/p1_paired_count_full"
```

Completed datasets have a validated result and a SHA-256 receipt. Missing, partial,
corrupted or job-only tasks are rerun. A batch contains at most five datasets;
completed datasets in a failed batch receive receipts when the parent is alive.
An abrupt power loss can require rerunning at most that batch's unreceipted work.
Workers use one BLAS thread. Changing worker count does not change seeds.
The OS writer lock prevents concurrent writes and is released on process death.
Never delete `run.lock` to force a second writer.

The full run writes `review.zip` with compact summaries, replicate metrics,
two PDF figures, warnings, source snapshots and provenance. Transfer that ZIP
for review; leave `tasks/` and `jobs/` on the server for resume verification.
There are no saved per-gene simulation matrices or millions of fitted models.

## What is the known truth?

The target is the difference between two coefficients for the SAME response,
not a genotype coefficient or a log-count coefficient. X0 = [1, binary exposure];
X1 adds three fixed quality covariates. One quality covariate has exposure
correlation 0.3 or 0.7; the other two are orthogonal to X0 and that covariate.
The full design has five columns. Sample sizes are 40, 80 and 160.

No quality covariate affects count generation. Exposure effects affect either
zero or 20% of genes, with balanced log2 effects of +/-0.75. EVERY gene is null
for the paired contrast, including genes with a genuine exposure effect.
Means span 2 to 2,000 expected counts; gene-specific NB dispersion is
0.02 + 0.3/sqrt(mean_count/10). These are specified benchmark parameters,
not a claim to be fitted human count parameters.

An optional `--parameters /path/to/count_generator_parameters.csv` uses the
existing `mean_count,dispersion` columns as generator inputs. Both must be
strictly positive and finite. A deterministic expression-stratified subset
preserves their pairing, and the file is copied and hashed in the new output.
This option does not reuse statistical tests or pretend to download new data.

**Null argument.** Within each binary exposure stratum, the joint distribution
of counts is invariant to sample permutations. Filtering by X0 and TMM/logCPM
transformation preserve this invariance. For each retained gene, the expected
transformed response is therefore constant within each exposure stratum and
belongs to the column space of X0. Any two full-rank WLS coefficient maps for
nested X0/X1 agree on that mean, so the paired mean contrast is exactly zero.
This remains true when the exposure changes a gene's mean. The weights depend
on the observed counts; the simulation tests the resulting inference rather
than assuming those estimated weights are independent or known.

Equal-depth scenarios have a fixed offset. In IID-depth scenarios, independent
lognormal depth multipliers (log SD 0.6, expected multiplier one) are redrawn
for EVERY dataset. Their null is marginal over these depth draws. They do not
claim calibration conditional on a particular fixed unequal-depth vector.
The TMM normalization factor and library total are estimated each time.
Before filtering/TMM, columns are randomly shuffled within exposure strata;
the original order is restored before voom and fitting. This symmetrizes any
reference-column tie handling, so the null argument does not depend on which
sample happened to appear first in a file. The permutation is seeded along
with the count draw. Sample metadata and model rows retain their fixed order.

The primary pipeline follows the original human analysis implementation:
`filterByExpr(design=X0)`, recompute filtered library totals, `calcNormFactors(method="TMM")`,
and `voom(adaptive.span=TRUE)`. See the official
[limma manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf)
for the transformation and precision-weight construction.

## Interpretation and report columns

- `fpr_050`: unadjusted rejection fraction among available paired tests, at 0.05.
- `fdp_planned_050`: mean FDP after BH over the planned gene family. Under this
  GLOBAL paired null it equals the probability of any rejection (FWER = FDR).
- `fdp_retained_050`: the corresponding BH result with the filtered family size.
  Unavailable retained tests occupy positions using p=1 for adjustment, but are
  never reported as valid nonsignificant tests.
- `coverage95`: coverage of zero among available intervals.
- `delivered_coverage95`: covered available intervals divided by all planned genes;
  interpret alongside filtering and availability, not as nominal interval coverage.
- FDR intervals are Wilson 95% binomial intervals across independent datasets.
  FPR/coverage Monte Carlo SEs use the dataset as the independent unit, not genes.
  At 1,000 datasets, the MC SE for a true FDR of 0.05 is about 0.0069.
- `pvalue_histogram.csv` is descriptive; correlated genes are not treated as
  independent observations in a formal uniformity test.

No success threshold is chosen after seeing results. Report the 0.05 target,
estimates and Monte Carlo intervals for ALL scenarios and BOTH weight modes.
Pipeline failures stop execution. An entirely unavailable arm is a failure,
not an FDR of zero. Warning messages are retained in `warnings.json`.

This experiment addresses end-to-end **global-null calibration**. It does not
estimate power for nonzero paired differences, strong FDR under mixed paired
nulls/alternatives, or calibration for arbitrary continuous exposure/nonlinear
quality effects. Such experiments need a separately justified transformed-response
truth. The specified count model does not include unmodeled donor/cage clustering.
An exposure effect in this benchmark is NOT a paired-contrast alternative.

## Tests

```bash
python -m unittest discover -s paired_count_calibration/tests -p 'test_*.py' -v
Rscript --vanilla paired_count_calibration/tests/test_kernel.R .
```

Use `--collect --output ACTUAL_DIRECTORY` to rebuild reports from verified receipts.
The experiment is separate from the existing manuscript tables and does not
overwrite the published-results snapshots.
