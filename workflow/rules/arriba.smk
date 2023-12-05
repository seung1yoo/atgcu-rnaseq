

rule star4arriba:
    input:
        fq1="anl/fastp/{sample}/{sample}_R1.fastq.gz",
        fq2="anl/fastp/{sample}/{sample}_R2.fastq.gz",
        idx=config["star_index"],
        gtf=config["genemodel_gtf"]
    output:
        bam="anl/arriba/{sample}/Aligned.out.bam",
        sortBam="anl/arriba/{sample}/Aligned.sortedByCoord.out.bam"
    log: "logs/star4arriba/{sample}.log",
    benchmark: "benchmarks/star4arriba/{sample}.tsv"
    message: "star4arriba is in progress on {wildcards.sample}"
    params:
        outprefix="anl/arriba/{sample}/"
    threads: 8
    conda: "envs/arriba.yml"
    shell:
        "STAR"
        " --runThreadN {threads}" 
        " --genomeDir {input.idx} --genomeLoad NoSharedMemory"
        " --readFilesIn {input.fq1} {input.fq2} --readFilesCommand zcat"
        " --outFileNamePrefix {params.outprefix}"
        " --outSAMtype BAM Unsorted SortedByCoordinate"
        " --outSAMunmapped Within --outBAMcompression 0"
        " --outFilterMultimapNmax 50 --peOverlapNbasesMin 10"
        " --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5"
        " --chimSegmentMin 10 --chimOutType WithinBAM HardClip"
        " --chimJunctionOverhangMin 10 --chimScoreDropMax 30"
        " --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1"
        " --chimSegmentReadGapMax 3 --chimMultimapNmax 50"
        " > {log} 2>&1"
        

rule arriba:
    input:
        bam=rules.star4arriba.output.bam,
        genome_fasta=config["genome_fasta"],
        gtf=config["genemodel_gtf"],
    output:
        fusions="anl/arriba/{sample}/fusions.tsv",
        discarded="anl/arriba/{sample}/fusions.discarded.tsv",
    log: "logs/arriba/{sample}.log"
    benchmark: "benchmarks/arriba/{sample}.tsv"
    message: "arriba is in progress on {wildcards.sample}"
    params:
        outprefix="anl/arriba/{sample}/"
    threads: 8
    conda: "envs/arriba.yml"
    shell:
        "arriba "
        " -x {input.bam}"
        " -o {output.fusions} -O {output.discarded}"
        " -a {input.genome_fasta} -g {input.gtf}"
        " -b {config[arriba_blacklist]}"
        " -k {config[arriba_known_fusions]}"
        " -t {config[arriba_known_fusions]}"
        " -p {config[arriba_protein_domains]}"
        " > {log} 2>&1"



