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
                          precomputed = decs,correct.all = T,
                          hvg.args = list(var.field="ratio",n=1000), # for use with cv2
                          PARAM=FastMnnParam(BSPARAM=BiocSingular::RandomParam()))

sce <- corrected$corrected

# make room
rm(corrected); rm(decs); gc()

sce$barcode <- colnames(sce)

colnames(sce) <- make.unique(colnames(sce))

# put old assay values back in ----------
add_uncorrected_expression <- function(assay = "logcounts", corrected=sce, raw_list=sces) {
  # get tbl of barcords/batches; will use for ensuring pairing of correct expression values by batch and barcode
  lkup <- colData(corrected) %>%
    as_tibble()
  
  # gets expression data from each rep
  all_logcounts <- map_df(raw_list, function(x) {x %>%
      assay(assay) %>%
      as_tibble(rownames = "gene_ID") %>%
      pivot_longer(-gene_ID, names_to = "barcode", values_to = "X") %>%
      filter(gene_ID %in% rownames(corrected))
  }, .id = "batch")
  
  # reshape into mat
  all_logcounts <- all_logcounts %>% pivot_wider(names_from = gene_ID, values_from = "X")
  
  # filter and reorder to match crrected
  all_logcounts <- left_join(lkup,all_logcounts)
  
  # sanity
  stopifnot(all(corrected$barcode == all_logcounts$barcode))
  stopifnot(all(corrected$batch == all_logcounts$batch))
  
  # corrected already has unique identifiers
  all_logcounts <- all_logcounts %>% mutate(identifier = colnames(corrected))
  
  # make a mat with identical features as the corrected mat
  mat <- dplyr::select(all_logcounts, -batch,-barcode) %>%
    column_to_rownames("identifier") %>%
    as.matrix() %>%
    .[,rownames(corrected)]
  
  # confirm they have identical rows and columns
  stopifnot(all(colnames(mat) == rownames(corrected)))
  stopifnot(all(colnames(corrected) == rownames(mat)))
  
  # add back in to the appropriate slot
  assay(corrected,assay) <- t(mat)  
  return(corrected)
}

# we don't really use the unspliced counts for anything yet, so won't add that back
sce <- add_uncorrected_expression(assay="logcounts")
sce <- add_uncorrected_expression(assay="counts")

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
