# Model Adjustment in Alzheimer Disease Transcriptomes

English installation and execution guide for version **1.2.2**. This pipeline runs on **Linux**.

[中文说明](README_CN.md)

## Overview

This repository evaluates how tissue-quality adjustment and cage modelling change expression coefficients, standard errors, and significant association sets. The main entry point is `run.py`.

The workflow includes public-data downloads, count aggregation, expression modelling, paired coefficient-difference tests, standard-error decomposition, dream convergence checks, simulations, calibration, and figures with CSV source data.

| Resource | Role |
|---|---|
| SEA-AD | Donor-level expression and tissue-quality adjustment in dorsolateral frontal cortex |
| GSE233208 | Animal-level spatial expression and alternative cage models |
| GSE118553 | Definition of the frontal candidate-gene panel |
| MGI | Human-to-mouse one-to-one ortholog mapping |

Analysis starts from publicly deposited count/UMI layers and processed microarray matrices. FASTQ alignment and raw microarray intensity processing are outside this workflow. Source URLs, eligibility rules, and processing details are documented in [Data sources](docs/DATA_SOURCES_EN.md) and implemented in `pipeline/` and `analysis/`. All eight detailed guides are available in English and Chinese; source-code comments remain in Chinese.

## Repository Layout

```text
paper1-from-public/
  run.py                   Main workflow
  run_dream_repair.py       Targeted review of dream warnings from an existing run
  config.json              Data release, eligibility rules, simulation settings, seeds
  environment.yml          Linux Conda environment
  requirements.txt         Python dependencies
  analysis/                R preparation, models, diagnostics, and simulations
  pipeline/                Python downloads, aggregation, calibration, and reporting
  tests/                   Regression and input-validation tests
  docs/                    Data, methods, output, and publication documentation
  CITATION.cff             Software citation metadata
  release_manifest.json    Checksums for the packaged release files
```

Raw data, prepared count matrices, installed R packages, and generated results are stored locally rather than in the repository.

## Choose a Run Mode

| Situation | Command option |
|---|---|
| First run, with data downloads allowed | `--task all` |
| Reuse an existing data directory; download only missing inputs if necessary | `--task all --data-work /path/to/existing/work` |
| Reanalyse complete prepared inputs without downloading data | `--task reanalyse --data-work /path/to/existing/work` |
| Install and verify R dependencies | `--task install-r` |

Run commands from the repository root, the directory containing `run.py`. Replace example paths with real Linux paths before execution.

## 1. First Run

### Create the environment

The supplied environment uses Python 3.12 and R 4.5.x. The R installer targets Bioconductor 3.22. Compilers and system libraries needed for R package installation are included in `environment.yml`.

```bash
cd /path/to/paper1-from-public
conda env create -f environment.yml
conda activate paper1-public
```

If the environment already exists, activate it instead of creating it again. An existing environment must provide the dependencies in `requirements.txt`, a compatible R installation, and the required build tools.

### Install R dependencies, then run the analysis

```bash
set -euo pipefail

python -u run.py \
  --task install-r \
  --work "$PWD/work_public" \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --library "$PWD/R_library" \
  2>&1 | tee install_r.log

python -u run.py \
  --task all \
  --work "$PWD/work_public" \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --library "$PWD/R_library" \
  --workers 1 \
  --resume \
  2>&1 | tee run.log
```

The shell stops if dependency installation or verification fails. The analysis itself does not automatically install packages. First-time installation and data acquisition require network access.

The full workflow checks dependencies and numerical regression tests before analysis. Downloaded files and prepared inputs remain in the work directory for subsequent use. Large source objects and intermediate count matrices require substantially more storage than the final review bundle; allow sufficient disk space and memory before starting.

## 2. Reuse Existing Data Without Downloading

Use the **complete work directory** from a compatible earlier run. A `review_bundle` is not a valid substitute: it deliberately excludes large count objects needed for analysis.

Set `DATA_WORK` to the existing data directory and `R_LIB` to the verified R library. Choose a new output directory for the new analysis.

```bash
cd /path/to/paper1-from-public
conda activate paper1
set -euo pipefail

DATA_WORK="/absolute/path/to/existing/work"
R_LIB="/absolute/path/to/existing/R_library"
NEW_WORK="$PWD/work_reanalysis_v1_2_2"

python -u run.py \
  --task reanalyse \
  --data-work "$DATA_WORK" \
  --work "$NEW_WORK" \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --library "$R_LIB" \
  --workers 1 \
  --resume \
  2>&1 | tee reanalysis.log
```

Here `paper1` is an example existing environment name; use the environment that contains your verified dependencies. This mode requires all four preparation stages and their recorded inputs to be complete. It checks input identities before refitting models, simulations, and calibration.

The source and destination work directories must be different and must not contain one another. Keep the original provenance records, checksums, and stage files with the existing data.

### Reuse data but allow missing files to be downloaded

In the command above, replace `--task reanalyse` with `--task all` and retain `--data-work`.

The workflow then:

1. Reuses prepared counts and metadata when their identities and provenance pass verification.
2. Rebuilds prepared inputs from verified raw caches when necessary.
3. Downloads only missing source files.

A file's presence alone is not sufficient for reuse. Missing integrity records or checksum mismatches require review; the pipeline does not silently accept an unverified cache.

## 3. Resume an Interrupted Run

Repeat the original command with the same code, configuration, paths, and `--resume` flag. Completed stages can be reused when their recorded signatures and outputs pass verification.

If code or analysis settings have changed, use a new output directory and point `--data-work` to the previous complete work directory. Do not delete stage records to bypass a signature mismatch.

Runtime depends on hardware, source-file access, R package compilation, and the number of model fits and simulation repeats. The full analysis is substantially longer than the preflight tests. Use the log and stage files to determine progress rather than expecting a fixed completion time.

## 4. Tests and Command-Line Options

Run the Python tests without loading the study datasets:

```bash
python -m unittest discover -s tests -v
```

Display all supported arguments:

```bash
python run.py --help
```

| Argument | Purpose |
|---|---|
| `--task` | Select the workflow or an individual task; default: `all` |
| `--work` | Output and working directory; default: `work_public` under the repository |
| `--data-work` | Existing work directory to reuse with `all` or `reanalyse` |
| `--reuse-work` | Alias for `--data-work` |
| `--config` | Configuration JSON; default: repository `config.json` |
| `--rscript` | Path to the Rscript executable |
| `--library` | R package library; default: repository `R_library` |
| `--workers` | Worker count; examples use one worker |
| `--resume` | Reuse validated completed stages |

Individual tasks also include `check`, `download`, `prepare`, `models`, `simulation`, `calibration`, `figures`, and `audit`. Downstream tasks require their upstream inputs. Use `all` or `reanalyse` for the complete analysis.

The standard preflight includes R numerical checks and small official-API dream tests. Passing these tests confirms that the tested components are usable; final completion is reported separately in `status.json`.

## 5. Find and Interpret the Results

All paths below are relative to the selected work directory.

| Location | Contents |
|---|---|
| `prepared/` | Aggregated counts, metadata, inclusion records, and design inputs |
| `results/` | Coefficients, SEs, p/q values, model diagnostics, and simulation results |
| `results/calibration/` | Paired-test calibration summaries and design geometry |
| `figures/` | Eight groups of analysis figures, with PDFs, 600-dpi PNGs, and source CSVs |
| `logs/` | Stage-specific process output |
| `stages/` | Stage records used for validation and resumption |
| `status.json` | Requested task status and audit summary |
| `review_bundle.zip` | Review package excluding large count objects |

Important result files include:

- `human_results.csv.gz`: human gene-context-contrast estimates in both model arms.
- `human_paired_changes.csv.gz`: paired coefficient differences and model-based/HC3 tests.
- `mouse_results.csv.gz`: mouse candidate-panel results and explicitly named testing families.
- `mouse_precision.csv.gz`: gene-level SE ratios and decomposition factors.
- `mouse_variance.csv.gz`: cage variance, conditional ICC, and fitting diagnostics.
- `mouse_model_status.csv`: availability and failures by specification and model.
- `mouse_dream_diagnostics.csv`: warnings requiring localization or review.
- `mouse_dream_retry_attempts.csv`: recorded optimizer retry results.
- `simulation_summary.csv`: bias, empirical SD, estimated SE, coverage, rejection, interval width, and availability.
- `calibration/calibration_summary.csv`: conditional calibration of the paired coefficient tests.

Files ending in `.csv.gz` are gzip-compressed CSV tables. For example, pandas can read them directly with `pd.read_csv("results/human_results.csv.gz")`; manual decompression is unnecessary.

### Exit codes

| Code | Meaning |
|---|---|
| `0` | The requested task completed |
| `2` | An input, dependency, execution, or validation check failed |
| `3` | The workflow finished but retained fit failures or unresolved warnings |

An `install-r` exit code of zero confirms dependency setup, not completion of the scientific analysis. Review the task named in `status.json` and the accompanying audit.

In the current simulation output, `failure_rate` records numerical failures, while `unavailable_rate` also includes structural non-estimability. A cage-fixed genotype effect cannot be estimated in a pure-genotype cage layout; those cases are recorded separately from numerical errors. Older manuscript source tables used different field names, so compare definitions before comparing columns.

## 6. Common Problems

| Message or symptom | What to check |
|---|---|
| Missing R packages | Run `--task install-r` with the same `--rscript` and `--library` intended for analysis |
| SSL or repository access error during installation | Resolve network/certificate access, then rerun dependency installation; do not interpret this as an analysis result |
| Prepared inputs missing | Confirm that `--data-work` points to a full work directory, not `review_bundle` |
| Checksum or provenance mismatch | Verify the source files and records; do not remove the validation records to force reuse |
| Code or parameter signature changed | Start a new `--work` directory and reuse verified data through `--data-work` |
| dream warning or exit code 3 | Inspect `mouse_dream_diagnostics.csv`, retry records, and `status.json`; see [Targeted dream review](docs/DREAM_REPAIR_EN.md) |
| More than one figure numbering scheme | Use the manuscript mapping document; code figure numbers do not directly replace manuscript figure numbers |

## 7. Relation to the Manuscript

The repository covers the core public-data analyses. Supplementary analyses of other cohorts and historical simulation designs have separate workflows and are not all rerun by this entry point.

Code figure numbers identify analysis modules, not the manuscript layout. Consult [Manuscript mapping](docs/PAPER_MAPPING_EN.md) before replacing or citing figures.

`results/manuscript_comparison.csv` retains an earlier manuscript's numerical reference values. In particular, its six mouse SE references predate the updated manuscript. Those constants do not enter model fitting or significance decisions. A mismatch in this historical comparison therefore does not by itself establish a mismatch with the current manuscript; compare the calculated result tables against the relevant manuscript version and testing family.

## Further Documentation

- [Data sources and processing](docs/DATA_SOURCES_EN.md)
- [Models, tests, and simulation definitions](docs/METHODS_EN.md)
- [Files and functions](docs/CODE_GUIDE_EN.md)
- [Output tables and fields](docs/OUTPUTS_EN.md)
- [Manuscript mapping and coverage](docs/PAPER_MAPPING_EN.md)
- [Targeted dream review](docs/DREAM_REPAIR_EN.md)
- [Release validation](docs/VALIDATION_EN.md)
- [GitHub upload guide](docs/GITHUB_UPLOAD_EN.md)

Use `CITATION.cff` for the available software citation metadata. Source datasets retain their providers' terms of use. The software license must be selected by the authors; public availability alone does not specify a reuse license.
