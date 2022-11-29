library(scater)
library(scran)
library(tidyverse)
library(batchelor)

# http://bioconductor.org/books/3.14/OSCA.multisample/integrating-datasets.html#quick-start

# -----------------------------
#  Get previously calculated HVG results
# -----------------------------

#dec_fls <- Sys.glob("results/downstream/single_sample_dimred/*.usa.filt.dimred.dec.rds")
dec_fls <- snakemake@input[["decs"]]

names(dec_fls) <- str_extract(dec_fls,"(?<=dimred\\/).+(?=\\.usa)")

decs <- dec_fls %>% map(read_rds)

# -----------------------------
#  Get the SCE opjects for each batch
# -----------------------------

#sce_fls <- Sys.glob("results/downstream/single_sample_clustering/*.usa.filt.dimred.clust.sce.rds")
sce_fls <- snakemake@input[["sces"]]

names(sce_fls) <- str_extract(sce_fls,"(?<=clustering\\/).+(?=\\.usa)")

sces <- sce_fls %>% map(read_rds)

# -----------------------------
#  Apply mnn via batchelor
# -----------------------------

# currently, correction is
corrected <- quickCorrect(sces,
                          precomputed = decs,
                          hvg.args = list(var.field="ratio",n=1000), # for use with cv2
                          PARAM=FastMnnParam(BSPARAM=BiocSingular::RandomParam()))

sce <- corrected$corrected

sce$barcode <- colnames(sce)

colnames(sce) <- make.unique(colnames(sce))

# ------------------- get cell info

cell_lookup <- read_tsv("data/cell_lookup.tsv.gz")

run_lookup <- read_tsv("data/gsm_run_lookup.tsv.gz") %>% dplyr::select(-Run) %>% distinct()

colData(sce) <- colData(sce) %>% 
  as_tibble(rownames = "unique.id") %>%
  left_join(run_lookup,by=c(batch = "Experiment")) %>%
  left_join(cell_lookup,by=c("barcode","orig.ident")) %>%
  column_to_rownames("unique.id") %>%
  DataFrame()

sce <- runTSNE(sce,dimred="corrected")

g_batches <- plotTSNE(sce,colour_by="batch") + facet_wrap(~colour_by)

#g_batches <- plotTSNE(sce,colour_by="batch")s

saveRDS(sce,snakemake@output[["sce"]])
saveRDS(g_batches,snakemake@output[["ggp"]])
write_tsv(g_batches$data,snakemake@output[["dat"]])
ggsave(snakemake@output[["png"]],g_batches)
