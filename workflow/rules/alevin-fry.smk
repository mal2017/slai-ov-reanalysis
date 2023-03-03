rule make_splici:
    """
    https://github.com/COMBINE-lab/alevin-fry
    https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/
    """
    input:
        #genome = refs_wf("results/repeatmasker/genome.fasta.masked"),
        genome = config.get("COMBINED_TE_GENOME_FA"),
        #gtf = rules.update_gtf_for_alevinfry.output.gtf,
        gtf = config.get("TXOME_AND_TES_GTF"),
    output:
        fasta = "results/transcriptome_splici_fl/transcriptome_splici_fl86.fa",
        t2g_3c = "results/transcriptome_splici_fl/transcriptome_splici_fl86_t2g_3col.tsv",
        t2g = "results/transcriptome_splici_fl/transcriptome_splici_fl86_t2g.tsv"
    resources:
        time=10,
        mem=4000,
        cpus=1
    params:
        read_length = 91,
        flank_trim_length = 5,
        odir="results/transcriptome_splici_fl/",
    conda:
        "../envs/eisar.yaml"
    priority:
        101
    script:
        "../scripts/make_splici.R"

rule make_splici_idx:
    input:
        fasta = rules.make_splici.output.fasta
    output:
        idx = directory("results/transcriptome_splici_idx")
    threads:
        24
    resources:
        time=480,
        mem=128000,
        cpus=24
    singularity:
        "docker://quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    priority:
        101
    shell:
        """
        salmon index -t {input.fasta} \
            -i {output.idx} -p {threads}
        """

rule download_or_copy2:
    output:
        temp(directory("results/fastqs/{sample}/"))
    resources:
        time=960,
        mem=20000,
        cpus=2
    threads:
        2
    #singularity:
    #    "docker://quay.io/biocontainers/sra-tools:3.0.0--pl5321hd0d85c6_1"
    priority:
        1
    conda:
        "../envs/sratools.yaml"
    shell:
        """
        mkdir -p {output}
        fasterq-dump --outdir {output}/ -e {threads} -p --split-files --include-technical {wildcards.sample}
        """

rule merge_by_experiment:
    input:
        lambda wc: expand("results/fastqs/{s}/",s=[x for x in SAMPLES if pep.get_sample(x).Experiment == wc.experiment])
    output:
        r1 = temp("results/merged_fastqs/{experiment}_1.fastq"),
        r2 = temp("results/merged_fastqs/{experiment}_2.fastq")
    params:
        r1 = lambda wc: expand("results/fastqs/{s}/{s}_1.fastq",s=[x for x in SAMPLES if pep.get_sample(x).Experiment == wc.experiment]),
        r2 = lambda wc: expand("results/fastqs/{s}/{s}_2.fastq",s=[x for x in SAMPLES if pep.get_sample(x).Experiment == wc.experiment])
    priority:
        99
    shell:
        """
        cat {params.r1} > {output.r1}
        cat {params.r2} > {output.r2}
        """

rule alevin_map:
    """
    https://alevin-fry.readthedocs.io/en/latest/getting_started.html#running-the-alevin-fry-pipeline
    library type choice: https://github.com/COMBINE-lab/salmon/discussions/674
    """
    input:
        idx = rules.make_splici_idx.output.idx,
        #fqdir =rules.download_or_copy2.output,
        r1 = rules.merge_by_experiment.output.r1,
        r2 = rules.merge_by_experiment.output.r2,
        tg = rules.make_splici.output.t2g
    output:
        dir = directory("results/alevin-fry/alevin/{experiment}")
    singularity:
        "docker://quay.io/biocontainers/salmon:1.5.2--h84f40af_0"
    threads:
        24
    resources:
        time=480,
        mem=128000,
        cpus=24
    params:
        #flags = lambda wc: pep.get_sample(wc.sample).ALEVIN_MAP_PARAM,
        flags = lambda wc: list(set([pep.get_sample(x).ALEVIN_MAP_PARAM[0] for x in SAMPLES if pep.get_sample(x).Experiment == wc.experiment]))[0], # this is a bit hacky - need to enforce one library type per run/biosample
        #r1 = "results/fastqs/{sample}/{sample}_1.fastq",
        #r2 = "results/fastqs/{sample}/{sample}_2.fastq",
    priority:
        100
    shell:
        """
        salmon alevin -i {input.idx} -p {threads} {params.flags} --sketch \
            -1 {input.r1} \
            -2 {input.r2} \
            --tgMap {input.tg} \
            -o {output.dir}
        """

rule alevin_fry_permit:
    input:
        map = rules.alevin_map.output.dir
    output:
        dir = directory("results/alevin-fry/quant/{experiment}")
    threads:
        1
    singularity:
        "docker://quay.io/biocontainers/alevin-fry:0.4.3--h7d875b9_0"
    resources:
        time=20,
        mem=4000,
        cpus=1
    priority:
        100
    shell:
        """
        alevin-fry generate-permit-list -d fw -k -i {input.map} -o {output.dir}
        """

rule alevin_fry_collate:
    input:
        pml = rules.alevin_fry_permit.output.dir,
        map = rules.alevin_map.output.dir
    output:
        touch("results/alevin-fry/quant/{experiment}.collate.done")
    threads:
        12
    singularity:
        "docker://quay.io/biocontainers/alevin-fry:0.4.3--h7d875b9_0"
    resources:
        time=20,
        mem=4000,
        cpus=12
    priority:
        100
    shell:
        """
        alevin-fry collate -t {threads} -i {input.pml} -r {input.map}
        """

rule alevin_fry_quant:
    input:
        tg3c = rules.make_splici.output.t2g_3c,
        quant = rules.alevin_fry_permit.output.dir,
        collation = rules.alevin_fry_collate.output,
    output:
        dir = directory("results/alevin-fry/quant_res/{experiment}")
    threads:
        12
    singularity:
        "docker://quay.io/biocontainers/alevin-fry:0.4.3--h7d875b9_0"
    resources:
        time=20,
        mem=4000,
        cpus=12
    priority:
        100
    shell:
        """
        alevin-fry quant -t {threads} \
            -i {input.quant} \
            -o {output.dir} \
            --tg-map {input.tg3c} \
            --resolution cr-like \
            --use-mtx
        """