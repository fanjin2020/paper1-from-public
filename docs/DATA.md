# Public Inputs

## Sources

| Resource | Use |
| --- | --- |
| GSE118553 / GPL10558 | Discovery panel and age sensitivity |
| SEA-AD DFC, release 2026-06-22 | Donor/subclass counts, clinical and quality models, broad-class mixtures |
| GSE233208 | Whole-section/cortex counts, animal identities and housing metadata |
| GSE48350 / GPL570 | Independent bulk expression sensitivity |
| GSE1297 / GPL96 | MMSE models and fixed signatures |
| GSE168137 | Mouse cortex estimated-count analyses |
| GSE157827 | All-recovered-nuclei donor aggregation |
| GSE92926 | Mouse microarray comparisons |
| GSE134151 | Deposited RPKM factorial analyses |
| GSE243292 | Author-labelled count reference and AD/PA comparisons |
| GSE167490 | Fresh target annotation and donor-level astrocyte analyses |
| GSE222007 | Fresh annotation, doublet sensitivity and cell-class interactions |
| MGI / HGNC | Orthology and explicit gene-identity mapping |

Executable URLs are in pipeline/sources.py, pipeline/seaad.py and pipeline/supplementary/. Public resource availability is external to this repository. Only authorized, publicly deposited objects are requested.

## Input Boundary

The pipeline starts from deposited raw count layers, estimated counts where explicitly declared, author-processed microarray expression/RPKM, and source metadata. It does not align FASTQ or regenerate arrays from CEL. Processed microarray/RPKM inputs are never relabelled as integer counts.

SEA-AD source release and donor inventory are explicit. Mouse cortex requires at least 100 eligible spots. Both main tissue summaries retain animal identities. Technical libraries are aggregated by their documented donor identifiers, not treated as independent donors.

External cohort analyses retain their assay-specific scales and covariates. Recovered-cell proportions are not tissue abundance estimates. The classifier reference is shared between designated analyses; this is not additional independent validation.

## Download Reuse

--raw-cache accepts directories containing source downloads. A requested object must match its registered URL and SHA256/size receipt. Filenames alone are insufficient. Conflicting versions or damaged receipts stop execution.

Verified immutable downloads may be hard-linked; otherwise they are copied. Downloads and their source receipts are the only cross-experiment data reuse. Derived RDS, donor NPZ, prepared matrices, inferred cell labels, result tables and simulation batches are not imported. Compressed-file extraction is performed in the current run.

Same-run donor or simulation checkpoints are allowed only with explicit continuation and matching experiment/code identities. Unfinished analysis files are not accepted merely because they exist.

## Storage

Keep output directories outside the code checkout. Source count objects and cell-level intermediates can be large; the final review export excludes them. No data deletion, remote publication or DOI registration is performed.

Software-library discovery reads installed package locations, not scientific results. The environment launcher never installs packages. Source-download SHA256 confirms local integrity against its receipt, not an independent publisher signature.
