# Code Guide

[English run guide](../README_EN.md) | [Chinese version](CODE_GUIDE.md)

## Main Entry Point

`run.py` orchestrates the complete workflow. `main` parses options, verifies Linux, obtains a directory lock, and schedules dependent stages. `fingerprint` records code and parameters; `stage` validates outputs before resuming; `rrun` starts an independent R process; `bundle` assembles the compact review package.

`prepare_human` aggregates SEA-AD and joins donor quality metadata. `prepare_mouse` reconstructs animal identity and reads Seurat counts. Verified input reuse is handled separately from preparation; neither preparation function substitutes old model results for count inputs.

## Python Modules

| File | Main functions | Responsibility |
|---|---|---|
| `pipeline/io.py` | `download`, `sha`, `dump`, `command` | HTTP files, integrity receipts, JSON, subprocess logs |
| `pipeline/sources.py` | `fetch_small` | Download enumerated source files to the current raw directory |
| `pipeline/sources.py` | `mouse_manifest` | Join GEO and supplementary metadata into animal identities and eligibility records |
| `pipeline/sources.py` | `unwrap_mouse` | Remove an extra gzip layer from the mouse R object |
| `pipeline/sources.py` | `quality_join` | Join official donor clinical and quality metadata |
| `pipeline/sources.py` | `read_geo_matrix`, `prepare_panel` | Parse GEO expression and sample fields into donor-region inputs |
| `pipeline/sources.py` | `map_orthologs` | Build one-to-one MGI mappings and record unavailable positions |
| `pipeline/seaad.py` | `sea_inventory`, `prepare_seaad` | Inventory the selected release and aggregate donor-subclass counts |
| `pipeline/seaad.py` | `h5_column`, `count_sum_sparse` | Decode HDF5 fields and sum sparse counts |
| `pipeline/seaad.py` | `collapse_genes`, `combine` | Merge duplicate symbols and combine donor matrices |
| `pipeline/calibration.py` | `geometry`, `evaluate` | Construct nested-design contrasts and evaluate differences, SEs, and F tests |
| `pipeline/calibration.py` | `mean_profiles`, `run` | Generate known-mean calibration conditions and performance summaries |
| `pipeline/figures.py` | `generate`, `save`, `ecdf` | Produce CSV source data, PDFs, and PNGs from computed results |
| `pipeline/reuse.py` | `collect`, `import_inputs` | Verify receipts/checksums and reuse count inputs without importing old fitted models |
| `pipeline/simulation_design.py` | `synthetic`, `matrices`, `prepare` | Generate exact allocations, design spaces, scenario grids, and seeds |
| `pipeline/audit.py` | `clinical_pairs`, `assess` | Check contrasts, denominators, failures, and historical manuscript anchors |
| `pipeline/cache.py` | See module documentation | Select verified prepared inputs, raw caches, or missing-file downloads |

Functions and classes retain their Chinese source documentation. SEA-AD preparation checkpoints contain aggregation work from the current downloaded inputs; the release does not ship precomputed count matrices.

## R Modules

| File | Responsibility |
|---|---|
| `analysis/preflight.R` | Install dependencies when requested; verify them in fresh processes; avoid an online Bioconductor version check when dependencies are already complete |
| `analysis/common.R` | Count validation, planned BH denominators, WLS differences, result extraction, and variance components |
| `analysis/diagnostics.R` | Weighted QR, SE implementation comparisons, dream interval checks, warning localization, optimizer agreement, and background-impact diagnostics |
| `analysis/panel.R` | Donor-level microarray models and frontal panel reconstruction |
| `analysis/prepare_mouse.R` | Seurat identity checks and whole-section/cortical count aggregation |
| `analysis/models.R` | Human arms and paired differences; five mouse models; precision and variance diagnostics |
| `analysis/simulation.R` | Gaussian four-model comparison on observed and synthetic layouts |
| `analysis/repair_dream.R` | Reconstruct affected inputs, verify default reproduction, retry optimizers, and moderate the merged background |

`common.R` is a shared library rather than a standalone workflow. R scripts receive their arguments from the Python entry points; opening an R script alone does not execute the complete pipeline.

## Targeted dream Entry Point

`run_dream_repair.py` is an optional migration/review entry point for complete v1.2.0 runs. It validates the earlier run and updates affected dream results in a separate output directory. New v1.2.2 analyses use `run.py`, which already contains the convergence checks. See [DREAM_REPAIR_EN.md](DREAM_REPAIR_EN.md).

## Tests

| File | Checks |
|---|---|
| `tests/test_calibration.py` | Independent least-squares coefficient maps, covariance, HC3, degeneracy, and known calibration means |
| `tests/test_inputs.py` | Identifiers, aggregation, GEO parsing, invalid gzip, cache integrity, partial HTTP responses, and documentation coverage |
| `tests/test_reanalysis.py` | Offline reuse, modified inputs, path containment, preparation versus modelling settings, duplicate endpoints, layout identities, and executed-code packaging |
| `tests/test_analysis.R` | Official EList versus explicit weighted-matrix interfaces, base-R WLS checks, invalid lists, zero weights, dimensions, and rank deficiency |
| `tests/test_cache.py` | Input reuse and cache selection |
| `tests/test_dream_repair.py` | Repair result integration |
| `tests/test_dream_retry.R` | Numerical acceptance rules and official-API server smoke tests |

Run Python tests without the study datasets:

```bash
python -m unittest discover -s tests -v
```

R regression tests run through `--task check` and the workflow preflight. The analysis tests require their declared R dependencies rather than silently skipping missing packages. Real-data analysis must wait until preflight passes.

## Making Changes

For a new model, document its estimand, input unit, testing family, and failure denominator alongside the implementation. Plotting code must display computed results rather than replace them with desired values.

Changing source URLs requires recording the source version and using a new work directory. Changes to cell annotations, age coding, or animal identity rules are analysis changes, not merely file-path fixes. Preserve the code and configuration associated with completed runs, and update release metadata when publishing modified code.
