library(AUCell)
library(scater)
library(scran)
library(tidyverse)
library(paletteer)

sce <- read_rds("~/amarel-matt/221123_grant_23_slaidina-pan-te-signature_01/results/downstream/integrated_by_tissue_and_sex_clustering/all.sce.rds")

sce <- sce[,!is.na(sce$active.ident)]

sigs <- read_rds("~/work/tetf_downstream/results/analysis/signatures/ourKD_gsea.tbl.rds") %>%
  filter(str_detect(comparison,"female_gonad")) %>%
  dplyr::select(kd,core_enrichment) %>%
  mutate(core_enrichment = map(core_enrichment,~str_split(.x,"/")[[1]])) %>%
  unnest(core_enrichment) %>%
  split(.,.$kd) %>%
  map(~filter(.x,core_enrichment %in% rownames(sce))) %>%
  map(pull,core_enrichment)

# --- signature
aggregated <- sumCountsAcrossFeatures(sce, sigs,
                                      exprs_values="reconstructed", average=T)
# -------- overview stuff

g_overview <- plotTSNE(sce,colour_by="active.ident")  + coord_equal() +   scale_color_paletteer_d("ggsci::default_igv")

g_ident <- plotTSNE(sce,colour_by="active.ident")  + facet_wrap(~colour_by) + coord_equal() +   scale_color_paletteer_d("ggsci::default_igv")

g_subs <- plotTSNE(sce,colour_by="data_subset")  + coord_equal() +   scale_color_paletteer_d("ggsci::default_igv")

g_batch <- plotTSNE(sce,colour_by="batch")  + coord_equal() +   scale_color_paletteer_d("ggsci::default_igv") + facet_wrap(~colour_by)


# ------- kd drivers

plotTSNE(sce, colour_by = "FBgn0000964",by_exprs_values = "reconstructed") + labs(caption = "tj expression (reconstructed after batch integration)") + 
  coord_equal()


# pan info ------------------------
plotTSNE(sce, colour_by = I(aggregated["pan",])) + labs(caption = "signature: pan-correlated TEs comprising core enrichment after knockdown") 

plotTSNE(sce, colour_by = "FBgn0085432",by_exprs_values = "reconstructed") + labs(caption = "pan expression (reconstructed after batch integration)") 

plotColData(sce, y=I(aggregated["pan",]), x="active.ident") + theme(axis.text.x = element_text(angle=45, hjust=1))
plotColData(sce,y = I(assay(sce)["FBgn0085432",]),  x="active.ident") + theme(axis.text.x = element_text(angle=45, hjust=1))
plotExpression(sce, features = "FBgn0085432", exprs_values = "reconstructed") + theme(axis.text.x = element_text(angle=45, hjust=1))




# ----- kds
plotTSNE(sce, colour_by = "FBgn0037698",by_exprs_values = "reconstructed") + labs(caption = "CG16779 expression (reconstructed after batch integration)")
plotTSNE(sce, colour_by = "FBgn0037698",by_exprs_values = "reconstructed") + labs(caption = "CG16779 expression (reconstructed after batch integration)")






plotColData(sce, y=I(aggregated["mamo",]), x="active.ident")
plotColData(sce, y=I(aggregated["ct",]), x="active.ident")

plotColData(sce, y=I(aggregated["CG16779",]), x="active.ident")



plotColData(sce, y=I(aggregated["awd",]), x="active.ident")





