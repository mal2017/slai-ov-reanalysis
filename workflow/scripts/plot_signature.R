library(AUCell)
library(scater)
library(scran)

sce <- read_rds("results/downstream/integrated_by_tissue_and_sex_clustering/all.sce.rds")


sigs <- read_rds("~/work/tetf_downstream/results/analysis/signatures/ourKD_gsea.tbl.rds") %>%
  filter(str_detect(comparison,"female_gonad")) %>%
  dplyr::select(kd,core_enrichment) %>%
  mutate(core_enrichment = map(core_enrichment,~str_split(.x,"/")[[1]])) %>%
  unnest(core_enrichment) %>%
  split(.,.$kd) %>%
  map(pull,core_enrichment)


#sigs <- rownames(sce)[!grepl("FBgn",rownames(sce))]

aggregated <- sumCountsAcrossFeatures(sce, sigs,
                                      exprs_values="reconstructed", average=TRUE)



sce$seurat_clusters <- as.factor(sce$seurat_clusters)



plotColData(sce, y=I(aggregated["pan",]), x="label")
plotTSNE()