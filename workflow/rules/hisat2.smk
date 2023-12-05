
if not Path(config["hisat2_index"]).exists():
    rule hisat2_index:
        input:
            fasta = config["genome_fasta"]
        output:
            directory(config["hisat2_index"])
        log: "logs/hisat2_index.log"
        benchmark: "benchmarks/hisat2_index.tsv"
        message: "hisat2_index is being run"
        threads: 1
        params:
            prefix = f"{config['hisat2_index']}/genome"
        wrapper:
            "master/bio/hisat2/index"

rule hisat2_align:
    input:
        reads=["anl/fastp/{sample}/{sample}_R1.fastq.gz",
               "anl/fastp/{sample}/{sample}_R2.fastq.gz"],
        idx=config["hisat2_index"]
    output:
        "anl/hisat2/{sample}/aligned.bam",
        "anl/hisat2/{sample}/{sample}.stats",
    log: "logs/hisat2_align/{sample}.log",
    benchmark: "benchmarks/hisat2_align/{sample}.tsv",
    message: "hisat2_align is in progress on {wildcards.sample}"
    params:
        extra="--summary-file anl/hisat2/{sample}/{sample}.stats",
    threads: 8
    wrapper:
        "master/bio/hisat2/align"

