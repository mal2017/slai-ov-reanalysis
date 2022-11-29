## Notes

Does not currently make use of subsamples, thus any line in the subsample table corresponds to a single
line in the sample table.

The index should be inspected before quantification and tmpfile deletion - depending on how the joint genome/TE fa and gtf were made it may exclude the TEs. This is often because the TE entries in the gtf lack "exon" tags or because the TEs are unintentionally duplicated in the gtf.


## Deps

`Snakemake`, `peppy`.

If not able to use the `--use-conda` flag, the following must be installed.

  - Seurat=4.3.0
  - tidyverse
  - SingleCellExperiment
  - zellkonverter
  - eisaR
  - celda
  - batchelor
  - PCAtools
  - fasterq-dump