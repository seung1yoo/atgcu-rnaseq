
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

#rule star_align:
#    input:
#        fq1="anl/fastp/{sample}/{sample}_R1.fastq.gz",
#        fq2="anl/fastp/{sample}/{sample}_R2.fastq.gz",
#        idx=config["star_index"]
#    output:
#        aln="anl/star/{sample}/Aligned.sortedByCoord.out.bam",
#        log="anl/star/{sample}/Log.out",
#        log_final="anl/star/{sample}/Log.final.out",
#        log_progress="anl/star/{sample}/Log.progress.out",
#        sj="anl/star/{sample}/SJ.out.tab",
#        reads_per_gene="anl/star/{sample}/ReadsPerGene.out.tab",
#        unmapped=["anl/star/{sample}/unmapped_R1.fastq.gz",
#                  "anl/star/{sample}/unmapped_R2.fastq.gz"],
#    log: "logs/star_align/{sample}.log",
#    benchmark: "benchmarks/star_align/{sample}.tsv",
#    message: "star_align is in progress on {wildcards.sample}"
#    params:
#        extra="--outSAMtype BAM SortedByCoordinate"
#        " --quantMode GeneCounts"
#        " --outSAMstrandField intronMotif"
#    threads: 8
#    wrapper:
#        "master/bio/star/align"

rule star_align:
    input:
        fastq_1="anl/fastp/{sample}/{sample}_R1.fastq.gz",
        fastq_2="anl/fastp/{sample}/{sample}_R2.fastq.gz",
        idx=config["star_index"]
    output:
        bam="anl/star/{sample}/Aligned.out.bam",
        sortBam="anl/star/{sample}/Aligned.sortedByCoord.out.bam",
        tranBam="anl/star/{sample}/Aligned.toTranscriptome.out.bam",
        readsPerGene="anl/star/{sample}/ReadsPerGene.out.tab",
        wigStr1="anl/star/{sample}/Signal.UniqueMultiple.str1.out.wig",
        wigStr2="anl/star/{sample}/Signal.UniqueMultiple.str2.out.wig",
        logFinal="anl/star/{sample}/Log.final.out",
        SJout="anl/star/{sample}/SJ.out.tab",
        unmapRead1="anl/star/{sample}/Unmapped.out.mate1",
        unmapRead2="anl/star/{sample}/Unmapped.out.mate2",
    wildcard_constraints:
        sample="[^/]+"
    params:
        outdir = 'anl/star/{sample}',
        prefix = 'anl/star/{sample}/',
    threads: 8
    message: "Running star on {wildcards.sample}"
    log: "logs/star_align/{sample}.log",
    benchmark: "benchmarks/star_align/{sample}.tsv",
    priority: 10
    conda: "../envs/star.yml"
    shell:
        "mkdir -p {params.outdir} && "
        "STAR --runMode alignReads"
        " --runThreadN {threads}"
        " --genomeDir {input.idx}"
        " --readFilesIn {input.fastq_1} {input.fastq_2}"
        " --readFilesCommand zcat"
        " --outFileNamePrefix {params.prefix}"
        " --outSAMtype BAM Unsorted SortedByCoordinate"
        " --outSAMstrandField intronMotif"
        " --quantMode TranscriptomeSAM GeneCounts"
        " --sjdbGTFfile {config[genemodel_gtf]}"
        " --outFilterType {config[star][outFilterType]}"
        " --outWigType wiggle"
        " --outWigStrand Stranded"
        " --outWigNorm RPM"
        " --outReadsUnmapped Fastx"
        #" --outFilterMultimapNmax {config[params_star_outFilterMultimapNmax]}"
        #" --alignSJoverhangMin {config[params_star_alignSJoverhangMin]}"
        #" --alignSJDBoverhangMin {config[params_star_alignSJDBoverhangMin]}"
        #" --outFilterMismatchNmax {config[params_star_outFilterMismatchNmax]}"
        #" --outFilterMismatchNoverLmax {config[params_star_outFilterMismatchNoverLmax]}"
        #" --alignIntronMin {config[params_star_alignIntronMin]}"
        #" --alignIntronMax {config[params_star_alignIntronMax]}"
        #" --alignMatesGapMax {config[params_star_alignMatesGapMax]}"
        #" --genomeFastaFiles {config[genomefasta]}"
        #" --outStd Log"
        #" --outSAMunmapped Within"
        #" > {output.bam}"

rule star_sortbam_idx:
    input:
        sortBam="anl/star/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        sortBamIdx="anl/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
    wildcard_constraints:
        sample="[^/]+"
    conda: "../envs/samtools.yml"
    shell:
        "samtools index {input.sortBam}"

rule parse_star_align_summary:
    input:
        expand("anl/star/{sample}/Log.final.out", sample=config["sample_names"])
    output:
        "anl/table/alignment_stats.star.xls"
    log: "logs/parse_star_align_summary.log"
    benchmark: "benchmarks/parse_star_align_summary.tsv"
    message: "parse_star_align_summary is being run"
    params:
        outprefix="anl/table/alignment_stats.star",
        target_dir="anl/star"
    shell:
        "python workflow/scripts/parser.py star"
        " --outprefix {params.outprefix}"
        " --target-dir {params.target_dir}"
        " > {log} 2>&1"

