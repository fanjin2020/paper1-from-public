# AdjScope-Tx Design Diagnostics

Version 0.2.0. The manuscript's design-geometry benchmark contains 48,000
Gaussian and 4,800 count datasets. It is separate from the core public-input
reproduction and its observed-layout count benchmark.

## Run

Copy `design_extension/`, `paper1_design_oneclick.py` and
`run_design_extension.sh` together into the server checkout. No GitHub clone,
package installation or new study-data download is performed.

```bash
conda activate paper1
bash run_design_extension.sh
```

The default is 36 independent workers, each with one R/BLAS thread. A new
`p1_design_MMDD_HHMMSS` sibling directory is printed at launch. Keep this exact
path for status and resume. Data preparation is serial. Wall time must be
estimated from the pilot on the target server, not from RAM capacity.

The launcher selects an existing Python with NumPy/SciPy; it does not enforce
the main reproduction's version pins. The selected versions and R library
paths are recorded and must remain unchanged for resume. R needs jsonlite,
edgeR, limma, variancePartition, BiocParallel, lme4, lmerTest and statmod.

Prepared inputs are discovered in sibling completed reproduction directories:
`core/prepared/mouse_units.csv`, `mouse_whole_counts.csv.gz`, and
`mouse_cortex_counts.csv.gz`. The newest eligible metadata file determines
automatic selection. Explicit selection is preferable when comparing runs:

```bash
bash run_design_extension.sh --data-work /path/to/completed/run/core \
  --r-library /path/to/existing/R_library
```

Repeat `--r-library` to combine existing libraries. These three inputs are
read-only and hashed. A preparation receipt is checked when present. Count
calibration is newly estimated; old fitted coefficients, p values, normalized
matrices, simulation outcomes and plotting extracts are not read.

## Profiles

| Profile | Layouts | Conditions | Count repetitions/cell | Count genes | Gaussian repetitions/cell |
|---|---:|---:|---:|---|---:|
| pilot | 12 | 4 | 3 | 600 expression-stratified | 100 |
| full | 12 | 8 | 50 | complete common background | 500 |

The pilot generates 144 count datasets and 4,800 Gaussian datasets. It tests
execution and exposes large effects; its repetitions cannot support the
selection rule. The full profile generates 4,800 count datasets and 48,000
Gaussian datasets. A count dataset is fitted by five models. Fifty repetitions
are an initial precision setting; inspect Monte Carlo uncertainty before
treating results as publication evidence. To use 200 repetitions, pass
`--count-reps 200` in a NEW full-profile output (19,200 count datasets).

```bash
bash run_design_extension.sh --profile full
```

Counts rerun filtering, TMM, voom/dream weight estimation, fitting, empirical
Bayes and BH using the existing five-model implementation. The count generator
uses empirical mean/dispersion/depth calibration from the whole-brain counts;
the eligible background is the intersection of the two mouse regions. The
genotype coefficient is a known library-offset log2 count-mean effect.

## Scientific Definitions

EWI is the within-cage exposure sum of squares divided by its total sum of
squares. It equals one minus the R-squared from cage indicators. This is a
standard information decomposition, not a newly discovered statistic.

Conditional EWI compares residual exposure information after baseline
covariates with information after adding cage indicators. Effective information
cages and the largest cage contribution describe concentration. No gene-level
estimated weights enter these pre-analysis measures.

In the simple model `Y = alpha + beta*x + gamma*xbar_cage + error`, with an
intercept, no other covariates and zero conditional mean error, the omitted
cage-mean bias of OLS is `gamma*(1-EWI)`. With additional covariates or a supplied
covariance matrix, the code computes the exact linear-estimator response to a
cage-mean shift. These are conditional Gaussian statements. They do not imply
an exact bias or FDR formula after nonlinear count transformations and weight
estimation. The Gaussian RE arm uses known covariance and is not dream.

An EWI of zero makes the within-cage genotype effect nonestimable. The tool
abstains instead of recommending FE or Mundlak at that boundary. An EWI of one
eliminates the specified genotype-related cage-mean shift in a balanced design;
it does not certify calibration for arbitrary errors or dependence.

## Prespecified Benchmark and Holdout

There are 72 synthetic animals per layout, with 12, 18 or 36 cages; the mixed
cage fraction is 0, 1/3, 2/3 or 1. Baseline covariates are age and sex. All
18-cage layouts are held out. The remaining layouts train the conditional rule.
This is held-out design generalization, not an independent biological cohort.

The conditions cross cage-mean amplitude 0/0.5 and latent cage SD 0/0.3. The
full profile adds the global null to the 10% nonnull condition. Twenty percent
of genes receive balanced signed cage-mean terms; nonnull genotype effects are
balanced +/-0.5. Full gene-level truth is regenerated deterministically rather
than exported as massive result tables.

The rule selects the highest training power among models passing a prespecified
empirical screen at the two nearest eligible training layouts. It requires at
least 20 usable training datasets, mean availability >=0.99, and a one-sided
95% Monte Carlo upper FDP-mean bound <=0.075. This tolerance is deliberately
reported separately from the BH level 0.05. Under the global null, every FDP is
an indicator and bounds are exact binomial Clopper-Pearson bounds. With nonnull
genes, bounds use a Student-t approximation for the mean dataset FDP and are
clipped to [0, 1]. When all FDP values are zero, a conservative zero-event
bound replaces the degenerate t interval. These bounds describe simulation uncertainty.

The rule is conditioned on the same generating cage-mean amplitude, cage SD
and nonnull fraction. These quantities are not known from real metadata.
Consequently, the real-data output contains design diagnostics, not a claimed
FDR forecast or universally valid model recommendation. Candidate selection
uses training metrics only; holdout outcomes evaluate the frozen rule. A new
run is required after changing the settings based on holdout findings.

## Monitor, Resume and Export

Replace the example output with the exact path printed at launch. Resume with
the same profile, input and scientific options; only worker count may change.

```bash
python paper1_design_oneclick.py --status --output /path/to/p1_design_run
bash run_design_extension.sh --resume --output /path/to/p1_design_run
python paper1_design_oneclick.py --collect --output /path/to/p1_design_run
```

Completed batches are reused only within this same new experiment after
SHA256/identity checks. Each count dataset has its own completion receipt.
Interrupted work without a valid receipt is rerun. Never delete `run.lock`.
Code, inputs or dependency changes require a fresh output directory. Stop the
foreground run with Ctrl+C; the Linux launcher terminates its worker groups.

Every finished run creates a short `p1_design_..._review.zip`. It contains
design summaries, analytic sensitivities, dataset-level performance, paired
cage-mean contrasts, a frozen training rule, holdout evaluation and HTML/Markdown
reports. No expression matrices, gene-by-gene fit tables or raw fit objects are
included in the review ZIP. Detailed diagnostic tables are retained for new
experiments; the compact report records their availability. Whole-arm failures
are prominently reported. A batch with zero usable inference
for an expected estimable arm stops the run early, so a broken model cannot
silently consume the entire simulation budget.

For metadata only, use `--design-csv /path/to/mouse_units.csv`. For a small
mechanism-only run without R or empirical input, use `--phase gaussian`.

## Testing and Provenance

```bash
python -m unittest discover -s design_extension/tests -v
Rscript --vanilla design_extension/tests/test_engine.R
```

The R engine retains the core fitted-model implementation and adds progress
records and diagnostic export. Unused shared-weight fixed-effect fits are
omitted in simulation-only execution. Source hashes are frozen in each run's
protocol. A completed legacy run can be summarized by `paper1_design_fix.py`
without refitting; resuming a legacy run requires its original source.

Not implemented: a universally valid EWI threshold, FDR Shapley attribution,
Hausman-triggered selection, an ML estimator, new biological cohorts, or an
installable R package. Unequal cage sizes and a systematic residual-dispersion
grid are not included in this first bounded study. Inference after data-dependent model choice would need
separate evaluation. The current deliverable is a Python/R command-line tool.
