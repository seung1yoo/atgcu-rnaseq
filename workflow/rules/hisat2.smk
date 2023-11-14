
if not Path(config["genome_hisat2_index"]).exists():
    rule hisat2_index:
        input:
            fasta = config["genome_fasta"]
        output:
            directory(config["genome_hisat2_index"])
        params:
            prefix = f"{config['genome_hisat2_index']}/genome"
            #prefix = config['genome_hisat2_index']
        log:
            "logs/hisat2_index.log"
        threads: 1
        wrapper:
            "master/bio/hisat2/index"

rule hisat2_align:
    input:
        reads=["anl/fastp/{sample}/{sample}_R1.fastq.gz",
               "anl/fastp/{sample}/{sample}_R2.fastq.gz"],
        idx=config["genome_hisat2_index"]
        #idx="src/hisat2_index"
    output:
        "anl/hisat2/{sample}/aligned.bam",
        "anl/hisat2/{sample}/{sample}.stats",
    log:
        "logs/hisat2/{sample}.log",
    params:
        extra="--summary-file anl/hisat2/{sample}/{sample}.stats",
    threads: 8
    wrapper:
        "master/bio/hisat2/align"

