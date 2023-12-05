
if not Path(config["star_index"]).exists():
    rule star_index:
        input:
            fasta=config["genome_fasta"],
            gtf=config["genemodel_gtf"]
        output:
            directory(config["star_index"]),
        log: "logs/star_index.log",
        benchmark: "benchmarks/star_index.log",
        message: "star_index is being run"
        threads: 1
        params:
            extra="",
            sjdbOverhang=config["star"]["sjdbOverhang"] # read_len - 1
        wrapper:
            "master/bio/star/index"

rule star_align:
    input:
        fq1="anl/fastp/{sample}/{sample}_R1.fastq.gz",
        fq2="anl/fastp/{sample}/{sample}_R2.fastq.gz",
        idx=config["star_index"]
    output:
        aln="anl/star/{sample}/Aligned.sortedByCoord.out.bam",
        log="anl/star/{sample}/Log.out",
        log_final="anl/star/{sample}/Log.final.out",
        log_progress="anl/star/{sample}/Log.progress.out",
        sj="anl/star/{sample}/SJ.out.tab",
        reads_per_gene="anl/star/{sample}/ReadsPerGene.out.tab",
        unmapped=["anl/star/{sample}/unmapped_R1.fastq.gz",
                  "anl/star/{sample}/unmapped_R2.fastq.gz"],
    log: "logs/star_align/{sample}.log",
    benchmark: "benchmarks/star_align/{sample}.tsv",
    message: "star_align is in progress on {wildcards.sample}"
    params:
        extra="--outSAMtype BAM SortedByCoordinate"
        " --quantMode GeneCounts"
        " --outSAMstrandField intronMotif"
    threads: 8
    wrapper:
        "master/bio/star/align"


