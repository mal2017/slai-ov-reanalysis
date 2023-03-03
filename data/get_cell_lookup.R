library(Seurat)
library(tidyverse)

x <- load("data/GSE162192_Seurat_objects.RData")


df <- tribble(~sra_sample, ~orig.ident, ~GSM,
        "adult1", "adult1", "GSM4946261",
        "adult2", "adult2","GSM4946262",
        "adult3", "adult4","GSM4946263",
        "adult4", "adult5","GSM4946264",
        "adult5", "adult6","GSM4946265")


get_cell_info <-  function(x) {
  x@meta.data %>% 
    as_tibble(rownames = "cell") %>% 
    mutate(barcode = str_remove(cell,"_.+"),suffix = str_extract(cell,"_\\d+")) %>%
    dplyr::select(cell,seurat_clusters,orig.ident,barcode,suffix) %>%
    mutate(active.ident = x@active.ident)
}

df_cells <- list(FC = FC_clusters, GC = GC_clusters, germarium = germarium_clusters) %>%
  map_df(get_cell_info, .id="data_subset")

df_runs <- read_csv("data/SraRunTable.txt") %>% dplyr::select(GSM=`GEO_Accession (exp)`,Run,Experiment) %>% distinct() %>% left_join(df)

write_tsv(df_cells,"data/cell_lookup.tsv.gz")
write_tsv(df_runs,"data/gsm_run_lookup.tsv.gz")
