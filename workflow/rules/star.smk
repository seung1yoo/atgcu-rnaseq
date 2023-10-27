
rule star_index:
    input:
        fasta=config["genome_fasta"],
        gtf=config["genemodel_gtf"]
    output:
        directory(config["genome_star_index"]),
    threads: 1
    params:
        extra="",
    log:
        "logs/star_index.log",
    wrapper:
        "master/bio/star/index"

rule star_align:
    input:
        fq1="anl/fastp/{sample}_R1.fastq.gz",
        fq2="anl/fastp/{sample}_R2.fastq.gz",
        idx=config["genome_star_index"]
    output:
        aln="anl/star/{sample}/aligned.bam",
        log="anl/star/{sample}/Log.out",
        log_final="anl/star/{sample}/Log.final.out",
        log_progress="anl/star/{sample}/Log.progress.out",
        sj="anl/star/{sample}/SJ.out.tab",
        #chim_junc="anl/star/{sample}/Chimeric.out.junction",
        reads_per_gene="anl/star/{sample}/ReadsPerGene.out.tab",
        unmapped=["anl/star/{sample}/unmapped_R1.fastq.gz",
                  "anl/star/{sample}/unmapped_R2.fastq.gz"],
    log:
        "logs/star/{sample}.log",
    params:
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts",
    threads: 8
    wrapper:
        "master/bio/star/align"

