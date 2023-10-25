
rule star_index:
    input:
        fasta=config["genome_fasta"],
        gtf=config["genemodel_gtf"]
    output:
        directory("src/star_index_genome"),
    threads: 1
    params:
        extra="",
    log:
        "logs/star_index_genome.log",
    wrapper:
        "master/bio/star/index"

rule star_align:
    input:
        fq1="anl/fastp/{sample}_R1.fastq.gz",
        fq2="anl/fastp/{sample}_R2.fastq.gz",
        idx="src/star_index_genome"
    output:
        aln="anl/star/{sample}/aligned.bam",
        log="anl/star/{sample}/Log.out",
        sj="anl/star/{sample}/SJ.out.tab",
        unmapped=["anl/star/{sample}/unmapped_R1.fastq.gz", "anl/star/{sample}/unmapped_R2.fastq.gz"],
    log:
        "logs/star/{sample}.log",
    params:
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate",
    threads: 8
    wrapper:
        "master/bio/star/align"
