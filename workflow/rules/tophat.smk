
def is_exists_tophat_index():
    prefix = config["genome_tophat_index"]
    for _ext in ["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.1.bt2"]:
        if not Path(f"{prefix}.{_ext}").exists():
            return 1
    return 0

if is_exists_tophat_index():
    rule tophat_index:
        input:
            fasta=config["genome_fasta"]
        output:
            multiext(config["genome_tophat_index"],
                     ".1.bt2",
                     ".2.bt2",
                     ".3.bt2",
                     ".4.bt2",
                     ".rev.1.bt2",
                     ".rev.2.bt2")
        log: "logs/tophat_index.log"
        benchmark: "benchmarks/tophat_index.tsv"
        params:
            extra=""
        threads: 8
        conda: "envs/tophat2.yml"
        shell:
            "bowtie2-build"
            " {params.extra}"
            " {input.fasta}"
            " {config[genome_tophat_index]}"

rule tophat_align:
    input:
        fq1="anl/fastp/{sample}/{sample}_R1.fastq.gz",
        fq2="anl/fastp/{sample}/{sample}_R2.fastq.gz",
        idx=multiext(config["genome_tophat_index"],
                     ".1.bt2",
                     ".2.bt2",
                     ".3.bt2",
                     ".4.bt2",
                     ".rev.1.bt2",
                     ".rev.2.bt2")
    output:
        bam="anl/tophat/{sample}/accepted_hits.bam",
        summay="anl/tophat/{sample}/align_summary.txt"
    log: "logs/tophat_align.{sample}.log"
    benchmark: "benchmarks/tophat_align.{sample}.tsv"
    params:
        idx=config["genome_tophat_index"],
        outdir="anl/tophat/{sample}",
        libraryType=lambda wildcards: config["samples"][wildcards.sample]["cuff_librarytype"]
    threads: 8
    conda: "envs/tophat2.yml"
    shell:
        "tophat --num-threads {threads}"
        " --output-dir {params.outdir}"
        " --library-type {params.libraryType}"
        " --rg-id {wildcards.sample}"
        " --rg-sample {wildcards.sample}"
        " --GTF {config[genemodel_gtf]}"
        " {params.idx} {input.fq1} {input.fq2}"
        

    

