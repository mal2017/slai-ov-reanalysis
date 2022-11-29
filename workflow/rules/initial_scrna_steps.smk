
rule make_r_objs:
    input:
        frydir = "results/alevin-fry/quant_res/{experiment}" #rules.alevin_fry_quant.output.dir
    output:
        sce = "results/downstream/raw_r_objs/{experiment}.usa.sce.rds",
    resources:
        time=5,
        mem=20000,
        cpus=2
    params:
        experiment = "{experiment}"
    conda:
        "../envs/make_r_objs.yaml"
    priority:
        3
    script:
        "../scripts/import_alevin_fry.R"

rule initial_filtering:
    """
    Applies low count outlier, mitochondrial, and doublet filtering per sample.
    """
    input:
        sce = "results/downstream/raw_r_objs/{experiment}.usa.sce.rds" #rules.make_r_objs.output.sce,
        #seur = rules.make_r_objs.output.seur,
    output:
        sce = "results/downstream/initial_filtering/{experiment}.usa.filt.sce.rds",
        g_doublet = "results/downstream/figs/gg/{experiment}_doublets.rds",
        g_miqc = "results/downstream/figs/gg/{experiment}_mito.rds",
        dat_doublet = "results/downstream/figs/tsv/{experiment}_doublets.tsv",
        dat_miqc = "results/downstream/figs/tsv/{experiment}_mito.tsv",
        png_doublet ="results/downstream/figs/png/{experiment}_doublets.png",
        png_miqc ="results/downstream/figs/png/{experiment}_mito.png",
    resources:
        time=240,
        mem=48000,
        cpus=2
    conda:
        "../envs/osca.yaml"
    priority:
        3
    script:
        "../scripts/initial_filtering.R"

rule sce_single_sample_dimred:
    """
    Finds HVGs, identifies number of PCs to use, and runs TSNE + UMAP per sample.
    """
    input:
        sce = rules.initial_filtering.output.sce
    output:
        sce = "results/downstream/single_sample_dimred/{experiment}.usa.filt.dimred.sce.rds",
        dec = "results/downstream/single_sample_dimred/{experiment}.usa.filt.dimred.dec.rds",
        g_hvg = "results/downstream/figs/gg/{experiment}_hvg.rds",
        g_elbow = "results/downstream/figs/gg/{experiment}_elbow.rds",
        dat_hvg = "results/downstream/figs/tsv/{experiment}_hvg.tsv",
        dat_elbow = "results/downstream/figs/tsv/{experiment}_elbow.tsv",
        png_hvg ="results/downstream/figs/png/{experiment}_hvg.png",
        png_elbow ="results/downstream/figs/png/{experiment}_elbow.png",
    resources:
        time=20,
        mem=20000,
        cpus=2
    conda:
        "../envs/osca.yaml"
    priority:
        3
    script:
        "../scripts/sce_single_sample_dimred.R"

rule sce_single_sample_clustering:
    input:
        sce = rules.sce_single_sample_dimred.output.sce
    output:
        sce = "results/downstream/single_sample_clustering/{experiment}.usa.filt.dimred.clust.sce.rds",
        g_sil = "results/downstream/figs/gg/{experiment}_sil.rds",
        dat_sil = "results/downstream/figs/tsv/{experiment}_sil.tsv",
        png_sil ="results/downstream/figs/png/{experiment}_sil.png",
        g_umap = "results/downstream/figs/gg/{experiment}_umap.rds",
        dat_umap = "results/downstream/figs/tsv/{experiment}_umap.tsv",
        png_umap ="results/downstream/figs/png/{experiment}_umap.png",
        g_tsne = "results/downstream/figs/gg/{experiment}_tsne.rds",
        dat_tsne = "results/downstream/figs/tsv/{experiment}_tsne.tsv",
        png_tsne ="results/downstream/figs/png/{experiment}_tsne.png",
    resources:
        time=20,
        mem=20000,
        cpus=2
    conda:
        "../envs/osca.yaml"
    priority:
        3
    script:
        "../scripts/sce_single_sample_clustering.R"


rule sce_integrate:
    input:
        sces = lambda wc: expand("results/downstream/single_sample_clustering/{s}.usa.filt.dimred.clust.sce.rds",s=EXPERIMENTS),
        decs = lambda wc: expand("results/downstream/single_sample_dimred/{s}.usa.filt.dimred.dec.rds",s=EXPERIMENTS)
    output:
        sce = "results/downstream/integrated_by_tissue_and_sex/all.sce.rds",
        ggp = "results/downstream/figs/gg/all_batches.rds",
        png =  "results/downstream/figs/png/all_batches.png",
        dat =  "results/downstream/figs/tsv/all_batches.tsv",
    resources:
        time=960,
        mem=256000,
        cpus=2
    conda:
        "../envs/osca.yaml"
    priority:
        4
    script:
        "../scripts/sce_integrate_by_tissue_and_sex.R"

rule sce_integrated_clustering:
    input:
        sce = rules.sce_integrate.output.sce
    output:
        sce = "results/downstream/integrated_by_tissue_and_sex_clustering/all.sce.rds",
        g_sil = "results/downstream/figs/gg/all_sil.rds",
        dat_sil = "results/downstream/figs/tsv/all_sil.tsv",
        png_sil ="results/downstream/figs/png/all_sil.png",
        g_tsne = "results/downstream/figs/gg/all_tsne.rds",
        dat_tsne = "results/downstream/figs/tsv/all_tsne.tsv",
        png_tsne ="results/downstream/figs/png/all_tsne.png",
    resources:
        time=60,
        mem=96000,
        cpus=2
    conda:
        "../envs/osca.yaml"
    priority:
        3
    script:
        "../scripts/sce_integrated_clustering.R"


# rule sce_integrated_dge:
#     input:
#         sces = lambda wc: expand("results/downstream/single_sample_clustering/{s}.usa.filt.dimred.clust.sce.rds",s=[x.sample_name for x in pep.samples if (x.tissue == wc.tissue) and (x.sex == wc.sex) and (x.sample_name in SAMPLES)]),
#         merged = rules.sce_integrated_clustering.output.sce,
#     output:
#         tsv = "results/downstream/integrated_by_tissue_and_sex_dge/{tissue}_{sex}.dge.tsv.gz",
#     resources:
#         time=60,
#         mem=64000,
#         cpus=2
#     conda:
#         "../envs/osca.yaml"
#     priority:
#         3
#     script:
#         "../scripts/sce_integrated_deg.R"
