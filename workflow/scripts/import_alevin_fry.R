library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
library(zellkonverter)

source("workflow/scripts/import_alevin_fry_func.R")

#proj <- "SRR13147913"
proj <- snakemake@params[["experiment"]]

#frydir <- "results/alevin-fry/quant_res/SRR13147913/"
frydir <- paste0(snakemake@input[["frydir"]])

x <- load_fry(frydir = frydir, which_counts = c("S","A")) # load spliced, which included ambig
u <- load_fry(frydir =frydir, which_counts = "U") # load unspliced

assay(x,"unspliced") <- counts(u) # combined into 1 sce

x@colData[["Run"]] <- proj

write_rds(x,snakemake@output[["sce"]])
#write_rds(sds,snakemake@output[["seur"]])
#zellkonverter::writeH5AD(sce = x,file = snakemake@output[["h5ad"]])
