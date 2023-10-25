
rule fastp_pe:
    #singularity: "docker://docker.io/staphb/fastp:0.23.4"
    input:
        sample=find_fastp_pe_input
        #fastq_1=lambda wildcards: config["samples"][wildcards.sample]["fastq_1"],
        #fastq_2=lambda wildcards: config["samples"][wildcards.sample]["fastq_2"]
    output:
        trimmed=["anl/fastp/{sample}_R1.fastq.gz", "anl/fastp/{sample}_R2.fastq.gz"],
        html="anl/fastp/{sample}.fastp.html",
        json="anl/fastp/{sample}.fastp.json"
    log:
        "logs/fastp/{sample}.log"
    threads: 10
    wrapper:
        "master/bio/fastp"

