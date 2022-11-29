library(eisaR)
library(Biostrings)
library(BSgenome)
library(stringr)
library(GenomicFeatures)

# see https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/

source("workflow/scripts/splici_func.R")

#gtf_path <- "results/tmp_af_txome.gtf"
gtf_path = file.path(snakemake@input[["gtf"]])
#genome_path <- "resources/plus-repeats.repeatmasked.fixednames.fasta"
genome_path = file.path(snakemake@input[["genome"]])

#read_length <- 91
read_length = snakemake@params[["read_length"]] # 91

#flank_trim_length <- 5
flank_trim_length = snakemake@params[["flank_trim_length"]] # 5

#output_dir = "tmp_test"
output_dir = snakemake@params[["odir"]]

dir.create(output_dir)

make_splici_txome(gtf_path=gtf_path, 
                  genome_path=genome_path, 
                  read_length=read_length, 
                  flank_trim_length=flank_trim_length, 
                  output_dir=output_dir)
