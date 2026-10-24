# Analysis Definitions

## Human Models

Eligible donor-context observations are aggregated from raw UMI counts. Base and quality-adjusted models use matched donors and gene backgrounds within context. Clinical status, conditional clinical status and ordinal Braak contrasts have distinct targets. The frontal panel is rediscovered, not imported or forced to a particular length. Missing or non-estimable positions remain explicit in planned testing families.

Paired coefficient differences use the covariance induced by fitting both models to the same observations. Model-based and HC3 covariance estimates are reported separately. Changing the BH family changes adjusted p values, not fitted coefficients or raw p values. Failure to reject a difference is not evidence of equivalence.

Matched sensitivity analysis crosses shared/native weights, ordinary/robust empirical Bayes and ordinal/categorical nuisance Braak coding. Categorical Braak changes the nuisance adjustment for the conditional clinical contrast; it does not redefine the ordinal Braak-slope target. Pairwise coefficient contrasts condition on estimated weights and are not empirical-Bayes-moderated paired tests.

## Mouse Models

Whole-section and cortex expression use the same eligible animals and a recomputed common gene background. Animal, cage and background counts are reported from the public inputs rather than embedded as substitutes for the data. Three covariate specifications per region give six settings. Seven arms comprise no cage, cage fixed effects, a common cage correlation, ordinary dream with shared/native weights and Mundlak dream with shared/native weights. Mundlak adds arithmetic cage means of baseline covariates. The target is the genotype coefficient.

The fixed-effect SE ratio is decomposed per gene into an unscaled-SE factor, posterior-SD factor and reweighting factor. Multiplication is checked at gene level; separate medians need not multiply to the median overall ratio. Weighted-centering identities are reported separately from the empirical arithmetic-mean Mundlak estimates.

Mixed-model boundary estimates, optimizer convergence and unavailable genes are recorded separately. Recoverable numerical fits are retried using the documented optimizer policy in `analysis/diagnostics.R`; successfully recovered raw rows are aligned before joint moderation. Unresolved fits are not converted to zero effects or p=1 without an explicit planned-family calculation.

## Count-to-Inference Validation

Eight scenarios cross latent cage log2 SD `{0, 0.3}`, genotype-related cage mean amplitude `{0, 0.5}` and nonzero-effect fraction `{0, 0.1}`. Each has 200 independent datasets with the observed animal allocation and the common whole-section background. Nonzero genotype effects are +/-0.5 log2 units. Dispersion and expression calibration are estimated from the prepared counts.

Every dataset repeats filtering, TMM normalization, voom weight estimation, fitting, empirical Bayes and BH for five model classes. This validation targets genotype coefficients, not a complete count-to-paired-difference procedure. Latent log2 cage SD is a generating parameter, not an observed count ICC.

FDR is the mean dataset-level false discovery proportion, including zero-discovery datasets. Power, coverage and availability retain their denominators in the exported tables. Entire-run failure is missing inference with explicit all-run FDR bounds. Planned-family and available-family BH are separate summaries. Monte Carlo uncertainty is calculated over independent datasets, not by treating genes as independent repetitions. All-null FDP/FWER intervals use the Wilson construction; non-null means use dataset-level Monte Carlo standard errors.

## Fixed-Design Benchmarks

The paired benchmark uses 5,000 responses per setting, three error laws, equal/fixed-gradient weights and varying sample size, partial correlation and leverage. Degenerate zero-difference maps are identified separately. Seeds are deterministic per design, error law and repetition, so batch partitioning does not change the paired streams.

For the nonzero-quality/zero-target-difference mean, let `c = d'Z`. Select the
lowest coordinate index within `256 * machine_epsilon * max(1, ||c||)` of the
smallest absolute coordinate of `c`, project that unit vector orthogonally to
`c`, and scale its resulting sample-level mean to RMS 0.5. This deterministic
tie rule prevents machine-roundoff differences between nominally zero coordinates
from selecting different simulated means. The coordinate, tolerance and projected
direction are recorded in `paired_geometry.csv`; the target-difference identity
is still checked for every batch. This changes neither the F-test formula nor
the paired covariance estimators.

The cage benchmark uses 5,000 repetitions per scenario over four allocations, two residual correlations, two cage-mean amplitudes and two target coefficients. No cage, fixed effects, common correlation, dream and Mundlak-dream share each generated dataset. Pure-genotype cages are structurally non-estimable for within-cage genotype coefficients in fixed/Mundlak models. The benchmark uses unmoderated inference; it is distinct from the count-to-inference experiment.

## Supplementary Public Analyses

Age sensitivity uses freshly prepared GSE118553 expression and identities, with overlap/quadratic designs and a baseline coefficient check. Quality variants use this run's SEA-AD counts, fixed-weight and design-specific voom alternatives, robust moderation and categorical nuisance Braak.

GSE48350 uses four explicitly distinguished age/scale variants. GSE1297 uses three MMSE contrasts, retaining a three-times-panel gene family and a separate six-signature family. Five external cohorts retain cohort-specific measurement scales, formulas and sample gates. Their primary focal family reserves 36 slots, including missing tests.

GSE167490 target labels are recomputed using the documented reference classifier, marker exclusions, confidence gates and donor-held-out checks. Technical replicates are collapsed by donor. The primary focal family has three slots; secondary models reserve twelve. Failed cell-class eligibility produces explicit unavailable rows, not missing family members. Recovered fractions remain conditional on sorting, QC and annotation.

GSE222007 annotations are regenerated from public counts. Scrublet scores are regenerated for every library with seed 290, with seeds 291/292 where the threshold diagnostic requires them. No stored cell labels or doublet scores enter a fresh run. Six cell-class contrasts share donors, matched count backgrounds and the specified second correlation estimate. Each count variant reserves six-times-panel tests; the combined exploratory family reserves twelve-times-panel tests. Annotation marker selection excludes endpoints and their mapped identities.

Seven broad SEA-AD classes use count-conserving aggregation. Complete matched donors supply the seven class models plus recovered-fraction and fixed-mean synthetic mixtures. The nine-times-panel gene family and seven recovered-composition tests remain separate.

Additional conditional Gaussian simulations use 1,737 responses per design/scenario and seed 319324. The separate three-model layout experiment uses 625 repetitions per scenario, 64 Gaussian background genes and scenario seeds derived from 3252026. Actual allocations come from current mouse metadata; synthetic allocations vary cage size, genotype balance and within-cage overlap. Infeasible layouts and non-estimable genotype coefficients remain explicit. These Gaussian experiments are not complete RNA-count simulations.

The supplementary Gaussian simulations preserve the row order in the freshly
prepared `core/prepared/mouse_units.csv`. Random animal errors are assigned by
this position, not by an animal-ID hash. `simulation_animal_order.csv` records
each draw position, animal identifier and covariates;
`simulation_input_order.json` records input/order hashes and seeds. Resuming a
partially completed simulation with a different order is rejected. An unchanged
seed alone does not identify an unchanged realized dataset when row order changes.
Neither recording nor validating this order sorts the data or changes the models.

## Reproducibility Controls

Seeds, a fresh-experiment identifier and scientific settings are part of the output protocol.
Only source download objects may be imported from another experiment. Prepared counts,
the panel, fitted models and simulation batches are recomputed in each new experiment.
The observed benchmark allocation is reconstructed from the newly prepared mouse identities.
Explicit --resume permits completed tasks from the same experiment only, after identity,
size and hash checks. Input and code identities are rechecked before final summaries.
Core stages also freeze the installed Python and R package versions. The top-level
runner checks the complete source identity before and after each subprocess.
All advertised whole-paper execution is subject to the explicit coverage gate.
The default integrated scope includes the central and supplementary analyses above,
20-section CSV exports and computed result displays. Coverage is a statement about
implemented workflows, not a numerical claim about a completed run or manuscript agreement.
Worker count is an execution setting, not a scientific parameter. Smoke tests reduce
repetitions and gene counts and are not manuscript evidence.

Matched task completion also checks numerical availability. An empirical model
with no usable genes or a wholly failed count-model run is not issued a successful
task receipt. Final completion requires an `ANALYSES_COMPLETE` summary. Individual
unavailable genes and prespecified structural non-estimability in the separate
cage benchmark remain explicit outcomes, not automatic errors.
