#!/usr/bin/env python
# vim: syntax=python tabstop=4 expandtab


rule run_fastqc:
    input:
        fastq_1=lambda wildcards: config["samples"][wildcards.sample]["fastq_1"],
        fastq_2=lambda wildcards: config["samples"][wildcards.sample]["fastq_2"]
    output:
        fastqc_1="analysis/fastqc/{sample}/{sample}_R1_fastqc.zip",
        fastqc_2="analysis/fastqc/{sample}/{sample}_R2_fastqc.zip"
    wildcard_constraints:
        sample="[^/]+"
    params:
        path="analysis/fastqc/{sample}/"
    threads: 4
    message: "Running fastqc on {wildcards.sample}"
    benchmark:
        "benchmarks/run_fastqc.{sample}.txt"
    shell:
        "fastqc -t {threads}"
        " --extract"
        " -o {params.path} {input.fastq_1} {input.fastq_2}"

