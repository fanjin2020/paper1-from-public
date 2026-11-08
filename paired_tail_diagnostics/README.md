# Paired-contrast tail diagnostics and weight controls

This is a separate, hypothesis-guided diagnostic experiment following the
completed 24-scenario count calibration. It does not change that code, its
receipts, or its results. It does not tune a new test to make the old FDR pass.

## One command

```bash
conda activate paper1
bash run_paired_tail_diagnostics.sh --profile full --workers 36 \
  --output "$HOME/adjscope/p1_paired_tail_full"
```

No downloads or installations. Use the existing manuscript R environment.
Full: six scenarios x 500 fresh datasets x 3,000 genes, six arms and two tests.
Pilot: six x 10 x 500. Smoke: six x one x 250. Override with `--reps`/`--genes`.
The full default gives 3,000 datasets, not another full 24,000-dataset run.
Its six arms mean runtime is not simply one eighth of the previous run.

```bash
bash run_paired_tail_diagnostics.sh --status --output "$HOME/adjscope/p1_paired_tail_full"
bash run_paired_tail_diagnostics.sh --resume --workers 36 --output "$HOME/adjscope/p1_paired_tail_full"
```

Each dataset has checksums for BOTH result.json and events.csv.gz. A job file,
partial file or damaged event table is not a completed dataset. Atomic writes,
an OS writer lock and deterministic seeds support resume. Output must not be
an existing count-calibration directory. The original paired kernel is sourced,
not replaced. Python uses its standard library; R uses edgeR, limma and jsonlite.

## Selected scenarios and known truth

Scenarios s01, s09, s11, s12, s15 and s23 cover n=40/80/160, equal/IID random
depths, exposure-quality correlations .3/.7, and zero/20% exposure-effect genes.
Selection was informed by the original results and is recorded in the protocol.
The main seed is 120261108, separated from the original simulation seeds.
Independent-weight counts use main seed + 50,000,000; Gaussian controls use
main seed + 100,000,000. The fixed covariate geometries and NB generating family
are retained. This is fresh simulation replication, not independent clinical validation.

All paired mean differences remain zero. Genuine exposure effects are not
alternatives to the paired-difference null. The new covariates have no count-mean
effect. The same within-exposure sample-order randomization and restoration is
used before normalization to symmetrize TMM reference ties. The IID-depth null
is marginal over the depths, not conditional on a fixed unequal-depth vector.
See the original `paired_count_calibration/README.md` for the null argument.

## Six arms

| Arm | Response | Base/full weights |
| --- | --- | --- |
| estimated_shared | Observed TMM logCPM | Same base-design voom weights |
| estimated_design | Same observed logCPM | Separate voom fits by model |
| independent_shared | Same observed logCPM | Base-design weights from independent NB counts |
| independent_design | Same observed logCPM | Separate weights from independent NB counts |
| fixed_unit | Same observed logCPM | All weights one |
| gaussian_unit | Independently generated Gaussian response | All weights one |

The auxiliary NB counts share the generator parameters, exposure effects and
latent depth multipliers with the observed data. Their counts and library totals
are independently realized. They use the observed retained gene set so gene
families and response matrices remain identical across the five count arms.
Thus independence is conditional on the latent depth and retained set: this
control removes outcome reuse in weight fitting, not all marginal dependence.
It is a simulation control, not a practical estimator requiring no extra data.

The fixed-unit arm removes precision weighting but retains TMM/log transformation
and count heteroskedasticity. It is not an oracle variance specification.
The Gaussian arm uses independent N(3 + .2*x, 1) responses on the same fixed
design and retained family. It checks HC3/reference-t behavior when the linear
Gaussian model and constant variance are correct. It is NOT compared to count
arms as if it shared their response. Its model-based statistic has the usual
Gaussian fixed-design t reference; HC3 remains a finite-sample sensitivity test.
Count attributes attached to Gaussian records are matched-gene labels, not
parameters generating that Gaussian response.

Each arm reports the unchanged unmoderated model-based and HC3 paired tests,
including t reference df=n-p_full. No empirical-Bayes variance moderation,
post-hoc df adjustment or threshold tuning is introduced.

## Outputs

Return only `review_tail.zip`; retain tasks/jobs on the server for resume.

- `summary.csv`: FPR, coverage, availability, planned/retained-family BH FDR
  with Wilson Monte Carlo intervals, plus Bonferroni family-wise error at .05/G.
  All twelve arm-test configurations remain visible.
- `tail_summary.csv`: strict counts/rates at 1e-6, 1e-5, .05/planned_genes,
  5e-5, 1e-4, 1e-3, .01 and .05. Rate/uniform-threshold ratios and MC SEs use
  independent datasets, not independent-gene assumptions. This is marginal
  calibration, not a formal independence test. Dataset-level probabilities of
  any exceedance have Wilson intervals. Zero observed tail events do not prove
  exact control; the empirical rate MC SE can then be zero while its uncertainty
  is not. Consult the any-event upper bound and increase replications if needed.
- `expression_strata.csv`: available genes and tail/rejection counts in generator
  mean bands <10, 10..<50, >=50. Empty strata retain zero counts and undefined rates.
  Means and dispersion are linked in this generator, so their effects cannot
  be causally separated by these strata alone.
- `paired_comparisons.csv`: within-dataset differences with MC SEs. Count-arm
  comparisons share the observed response; HC3/model-based comparisons share
  their respective arm. Gaussian-vs-count paired comparisons are excluded.
- `tail_events.csv.gz`: extreme/false-discovery gene records with both fitted
  coefficients, delta, both SEs, t/df/p, both q values, maximum leverage, weight
  ratio, generator mean/dispersion, observed count mean and zero fraction.
  Records include all BH discoveries, p<=1e-4, and the minimum-p gene per arm/test.
  These are selected diagnostic records, NOT a full test family for redoing BH.
  A safety cap of 200 records per arm/test/dataset retains the smallest p values;
  `events_truncated` explicitly reports excluded records. Aggregate counts
  always use every tested gene, even if a diagnostic export is capped.
- `replicate_metrics.csv.gz`: compact per-dataset metrics for independent review.
- `protocol.json`, `environment.json`, `warnings.json`, source and checksum manifest.

New histogram counts use `fuzz=0`; strict threshold counts do not depend on
histogram bins. Old histograms and production results are not rewritten.

## Interpretation

If Gaussian/HC3 also has inflated extreme tails, estimated voom weights are
not necessary for that behavior. If independently estimated weights reduce
inflation, outcome reuse in weight estimation is implicated under the tested
depth/filtering conditions. Neither comparison alone proves a unique cause.
If marginal tails look controlled but family error remains elevated, investigate
joint dependence rather than equating pointwise coverage with BH control.
Bonferroni needs no independence assumption: convincing inflation of its
global-null family-wise error cannot be explained solely by BH's dependence
conditions. A noninflated Bonferroni result alone does not prove every marginal
p threshold valid or identify a unique explanation for BH inflation.

No arm is automatically recommended after seeing the same benchmark results.
These controls diagnose the failure; a replacement inferential method needs
its own justified null calibration and held-out evaluation. Nonzero paired
effect power, strong FDR under mixed paired nulls, and clinical replication
are not added by this module.
