# Files excluded from the initial GitHub upload

The repository includes all files from `input/` and all `output/` files smaller than GitHub's 100 MB per-file limit. The following regenerable or oversized data are intentionally excluded:

- `intermediate/` (about 1.5 GB): regenerable workflow state.
- `renv/library/` and other local renv caches (about 1.3 GB): restore from `renv.lock`.
- `sc_reference_preprocessing/Wu et al.BRCA/wu2021_brca_seurat.rds` (about 798 MB).
- `output/CRC1/02_morphology_cluster/CRC1_raw_SME_normalizeA.mtx` (about 820 MB).
- Eleven CRC1 inferCNV object checkpoints between about 121 MB and 364 MB each.

These files should be distributed through a data repository or Git LFS if they need to be versioned later.