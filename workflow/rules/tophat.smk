
def is_exists_tophat_index():
    prefix = config["tophat_index"]
    for _ext in ["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.1.bt2"]:
        if not Path(f"{prefix}.{_ext}").exists():
            return 1
    return 0

if is_exists_tophat_index():
    rule tophat_index:
        input:
            fasta=config["genome_fasta"]
        output:
            multiext(config["tophat_index"],
                     ".1.bt2",
                     ".2.bt2",
                     ".3.bt2",
                     ".4.bt2",
                     ".rev.1.bt2",
                     ".rev.2.bt2")
        log: "logs/tophat_index.log"
        benchmark: "benchmarks/tophat_index.tsv"
        message: "tophat_index is being run"
        params:
            extra=""
        threads: 1
        conda: "../envs/tophat.yml"
        shell:
            "bowtie2-build"
            " {params.extra}"
            " {input.fasta}"
            " {config[tophat_index]}"
            " > {log} 2>&1"

rule tophat_align:
    input:
        fq1="anl/fastp/{sample}/{sample}_R1.fastq.gz",
        fq2="anl/fastp/{sample}/{sample}_R2.fastq.gz",
        idx=multiext(config["tophat_index"],
                     ".1.bt2",
                     ".2.bt2",
                     ".3.bt2",
                     ".4.bt2",
                     ".rev.1.bt2",
                     ".rev.2.bt2")
    output:
        bam="anl/tophat/{sample}/accepted_hits.bam",
        summay="anl/tophat/{sample}/align_summary.txt"
    log: "logs/tophat_align/{sample}.log"
    benchmark: "benchmarks/tophat_align/{sample}.tsv"
    message: "tophat_align is in progress on {wildcards.sample}"
    params:
        idx=config["tophat_index"],
        outdir="anl/tophat/{sample}",
        libraryType=lambda wildcards: config["samples"][wildcards.sample]["cuff_librarytype"]
    threads: 8
    conda: "../envs/tophat.yml"
    shell:
        "tophat --num-threads {threads}"
        " --output-dir {params.outdir}"
        " --library-type {params.libraryType}"
        " --rg-id {wildcards.sample}"
        " --rg-sample {wildcards.sample}"
        " --GTF {config[genemodel_gtf]}"
        " {params.idx} {input.fq1} {input.fq2}"
        " > {log} 2>&1"
        
rule parse_tophat_align_summary:
    input:
        expand("anl/tophat/{sample}/align_summary.txt", sample=config["sample_names"])
    output:
        "anl/table/alignment_stats.tophat.xls"
    log: "logs/parse_tophat_align_summary.log"
    benchmark: "benchmarks/parse_tophat_align_summary.tsv"
    message: "parse_tophat_align_summary is being run"
    params:
        outprefix="anl/table/alignment_stats.tophat",
        target_dir="anl/tophat",
        samples=expand("{sample}", sample=config['sample_names'])
    shell:
        "python workflow/scripts/parser.py tophat"
        " --outprefix {params.outprefix}"
        " --target-dir {params.target_dir}"
        " --sample-id-s {params.samples}"
        " > {log} 2>&1"


