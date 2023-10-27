
rule hisat2_index:
    input:
        fasta = config["genome_fasta"]
    output:
        directory(config["genome_hisat2_index"])
    params:
        prefix = config["genome_hisat2_index"]
    log:
        "logs/hisat2_index.log"
    threads: 1
    wrapper:
        "master/bio/hisat2/index"

rule hisat2_align:
    input:
        reads=["anl/fastp/{sample}_R1.fastq.gz",
               "anl/fastp/{sample}_R2.fastq.gz"],
        idx=config["genome_hisat2_index"]
    output:
        "anl/hisat2/{sample}/aligned.bam",
    log:
        "logs/hisat2/{sample}.log",
    params:
        extra="",
    threads: 8
    wrapper:
        "master/bio/hisat2/align"
