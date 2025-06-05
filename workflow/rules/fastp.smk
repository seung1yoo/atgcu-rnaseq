
rule fastp_pe:
    #singularity: "docker://docker.io/staphb/fastp:0.23.4"
    input:
        sample=find_fastp_pe_input
        #fastq_1=lambda wildcards: config["samples"][wildcards.sample]["fastq_1"],
        #fastq_2=lambda wildcards: config["samples"][wildcards.sample]["fastq_2"]
    output:
        trimmed=["anl/fastp/{sample}/{sample}_R1.fastq.gz", 
                 "anl/fastp/{sample}/{sample}_R2.fastq.gz"],
        html="anl/fastp/{sample}/{sample}.fastp.html",
        json="anl/fastp/{sample}/{sample}.fastp.json"
    log: "logs/fastp_pe/{sample}.log"
    benchmark: "benchmarks/fastp_pe/{sample}.tsv"
    message: "fastp_pe is in progress on {wildcards.sample}"
    threads: 10
    wrapper:
        "master/bio/fastp"

rule parse_fastp_raw:
    input:
        expand("anl/fastp/{sample}/{sample}.fastp.json", sample=config["sample_names"])
    output:
        "anl/table/sequencing_stats.raw.xls"
    log: "logs/parse_fastp_raw.log"
    benchmark: "benchmarks/parse_fastp_raw.tsv"
    message: "parse_fastp_raw is being run"
    params:
        outprefix="anl/table/sequencing_stats.raw",
        target_dir="anl/fastp",
        target_tag="before_filtering",
        samples=expand("{sample}", sample=config['sample_names'])
    shell:
        "python workflow/scripts/parser.py fastp"
        " --outprefix {params.outprefix}"
        " --target-dir {params.target_dir}"
        " --target-tag {params.target_tag}"
        " --sample-id-s {params.samples}"
        " > {log} 2>&1"

rule parse_fastp_clean:
    input:
        expand("anl/fastp/{sample}/{sample}.fastp.json", sample=config["sample_names"])
    output:
        "anl/table/sequencing_stats.clean.xls"
    log: "logs/parse_fastp_clean.log"
    benchmark: "benchmarks/parse_fastp_clean.tsv"
    message: "parse_fastp_clean is being run"
    params:
        outprefix="anl/table/sequencing_stats.clean",
        target_dir="anl/fastp",
        target_tag="after_filtering",
        samples=expand("{sample}", sample=config['sample_names'])
    shell:
        "python workflow/scripts/parser.py fastp"
        " --outprefix {params.outprefix}"
        " --target-dir {params.target_dir}"
        " --target-tag {params.target_tag}"
        " --sample-id-s {params.samples}"
        " > {log} 2>&1"

