# Re-export a Completed Experiment

This operation reads one completed public-input experiment. It does not refit models,
repeat simulations, download data or install software. It is distinct from a fresh
reproduction, which must still recompute all scientific analyses.

It now generates `reporting/Supplementary_Tables.xlsx` with 20 worksheets and 112
blocks, including the seven detailed-result blocks absent from the archived
v2.1.0 workbook map. Both completed legacy 105-block runs and verified portable
collections can be inputs. Original scientific results and completion records
remain unchanged. The exporter version is `publication_exports_20261026_v3`.

```bash
conda activate paper1
cd /home/research/git/repository/paper1-from-publicV1
python paper1_reexport_results.py \
  --run /home/research/git/repository/paper1_complete_from_public_20261022_v2
```

The script selects the interpreter recorded in `software.json` and checks only the
four export dependencies: NumPy, pandas, SciPy and Matplotlib. It does not require
openpyxl, an R installation, or the single-cell preparation packages. An unavailable
recorded interpreter can be replaced explicitly with `--python /existing/python`
only when its reporting versions match the record.

The output is a new sibling directory named
`p1_review_YYYYMMDD_HHMMSS`, for example `p1_review_20261024_153012`.
The collector uses the same short naming convention. An explicit `--output` overrides
the default; an existing directory is never overwritten. On Windows, transfer to a short
parent such as `C:\p1\` and do not nest the result folder inside another copy of its name.
The last terminal line beginning `REVIEW_DIRECTORY` gives its exact path.
Transfer that entire directory after `COLLECTION_COMPLETE.json` appears. No ZIP
is created. Do not send only the figures or only the reporting directory: the
package also carries detailed fits, count examples, unavailable-test records and
the source/environment receipts required to check the manuscript.

```bash
python paper1_reexport_results.py --verify /path/to/the/exported/directory
```

## Definitions

- ICC tables join dream SE to baseline SE by gene, specification and weight policy.
  The bins are <0.01, 0.01-0.10, 0.10-0.20 and >=0.20, with left-closed boundaries.
- Native fixed/base precision has unscaled-SE and posterior-SD factors. Shared-weight
  precision has separate design, posterior-SD and reweighting factors. They are not interchangeable.
- BH exports include the human family and separate two- and three-arm mouse families,
  with planned and available-only denominators kept in different columns.
- Main Figure 4 has four panels. Supplementary Figures S1-S12 follow the manuscript's
  quality, astrocyte, consistency, ICC, rejection, continuous-change, design, calibration,
  workflow, selected-change, donor-calibration and cancellation/power ordering.
- Main Figure 2 labels planned positions per family. Figure 5 marks the n=80,
  ordinary-leverage, uniform-weight design with partial correlation 0.3 and labels
  panels A-F. Figure 6 labels rows A-C and annotates Mundlak FDR estimates near 5%.
- `block_sources.csv` has one explicit source per block. Source tables can contain
  additional diagnostic fields. S2.3 checks within-run fitted numerical identities.
  S16.3 documents precision-component definitions.

- S1.2 links to public resource pages and input definitions, including SEA-AD DFC
  release 2026-06-22. Raw download receipts remain in the run archive.
- S6.6 and S10.2 provide context-specific sex counts and CAMERA results.
  S20.13-S20.16 provide full count-pipeline details, and S20.17 expands every
  simulation job from the recorded protocol without changing its definition.
- Parameter-table input filenames identify generated files. Software locations
  are displayed by package name. The CSV/JSON sources retain original paths.
  The XLSX writer uses the Python standard library and adds no package downloads.

To rebuild only the workbook from an already generated 112-block reporting
directory, run `python -m pipeline.workbook --reporting /path/to/reporting`.
For a 105-block directory, use the full re-export command above to first generate
the missing blocks from the completed experiment. Workbook generation requires
the complete indexed source tables and verifies their hashes and row counts.

- S4.2 summarizes all finite paired clinical log components. `n_decomposable` is
  the distribution denominator; `n_selected` counts base q < 0.05 separately.
- S9 leaves disease-direction comparisons missing for interaction coefficients such
  as `case:treat`, unavailable estimates and unknown discovery directions. The
  `direction_comparison_status` column explains why; testability, beta, p and q stay unchanged.
- S19.4 uses actual result layouts for `simulated_five_model` and
  `simulated_three_model`. `simulated_here` is an alias of the five-model flag,
  not of feasibility.
- S15.1 `coefficient_join_difference` is absolute. The separate
  `coefficient_join_signed_difference` records dream beta minus variance-fit beta.

## Provenance and Limits

The original experiment is read-only, including its code identity and completion records.
The exporter has its own code identity in `export_provenance.json`; its source is under
`export_source_snapshot`. It is not mislabeled as the source that fitted the experiment.
Raw inputs and large temporary simulation objects are not copied.
For corrected publication definitions, use `reporting/block_sources.csv` and its
mapped tables. Immutable scientific-source copies can retain earlier derived labels;
their preservation is for provenance, not an alternative publication definition.

Keep all three S18 blocks from the same completed experiment. Their Monte Carlo
realizations depend on sample order as well as seeds, scenario settings and
software. Re-exporting preserves the fitted values and does not regenerate random
draws or replace them with publication targets. Successful export means the files
and defined summaries were generated; manuscript values must be checked against
that specific experiment.
