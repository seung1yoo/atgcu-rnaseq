
rule cufflinks:
    input:
        bam=rules.tophat_align.output.bam
        #bam=rules.star_align.output.sortBam
    output:
        fpkmGene="anl/cufflinks/{sample}/genes.fpkm_tracking",
        fpkmTran="anl/cufflinks/{sample}/isoforms.fpkm_tracking"
    params:
        outdir='anl/cufflinks/{sample}',
        libraryType=lambda wildcards: config["samples"][wildcards.sample]["cuff_librarytype"]
    log: "logs/cufflinks/{sample}.log"
    benchmark: "benchmarks/cufflinks/{sample}.tsv"
    message: "cufflinks is in progress on {wildcards.sample}"
    threads: 8
    conda: "../envs/cufflinks.yml"
    shell:
        "cufflinks -p {threads}"
        " --multi-read-correct"
        " --frag-bias-correct {config[genome_fasta]}"
        " --library-type {params.libraryType}"
        " -G {config[genemodel_gtf]}"
        " -o {params.outdir}"
        " {input.bam}"
        " > {log} 2>&1"

rule parse_cufflinks:
    input:
        expand("anl/cufflinks/{sample}/genes.fpkm_tracking", sample=config["sample_names"])
    output:
        "anl/table/abundance.cufflinks.genes.fpkm.all.xls",
        "anl/table/abundance.cufflinks.genes.fpkm.exp.xls",
        "anl/table/abundance.cufflinks.genes.stats.exp.xls",
        "anl/table/abundance.cufflinks.genes.stats.coexp.xls",
        "anl/table/abundance.cufflinks.isoforms.fpkm.all.xls",
        "anl/table/abundance.cufflinks.isoforms.fpkm.exp.xls",
        "anl/table/abundance.cufflinks.isoforms.stats.exp.xls",
        "anl/table/abundance.cufflinks.isoforms.stats.coexp.xls"
    log: "logs/parse_cufflinks.log"
    shell:
        "python workflow/scripts/parser.py cufflinks"
        " --outprefix anl/table/abundance.cufflinks"
        " --gene2desc {config[anno_gene2desc]}"
        " --tran2desc {config[anno_tran2desc]}"

rule cuffquant:
    input:
        #bam=rules.tophat_align.output.bam
        bam=rules.star_align.output.sortBam
    output:
        cxb="anl/cuffquant/{sample}/abundances.cxb"
    log: "logs/cuffquant/{sample}.log"
    benchmark: "benchmarks/cuffquant/{sample}.tsv"
    message: "cuffquant is in progress on {wildcards.sample}"
    params:
        outdir='anl/cuffquant/{sample}',
        libraryType=lambda wildcards: config["samples"][wildcards.sample]["cuff_librarytype"]
    threads: 8
    conda: "../envs/cufflinks.yml"
    shell:
        "cuffquant -p {threads}"
        " -o {params.outdir}"
        " --multi-read-correct"
        " --frag-bias-correct {config[genome_fasta]}"
        " --library-type {params.libraryType}"
        " {config[genemodel_gtf]}"
        " {input.bam}"
        " > {log} 2>&1"

rule cuffnorm:
    input:
        expand("anl/cuffquant/{sample}/abundances.cxb", sample=config["sample_names"])
    output:
        "anl/cuffnorm/genes.fpkm_table",
        "anl/cuffnorm/isoforms.fpkm_table"
    log: "logs/cuffnorm.log"
    benchmark: "benchmarks/cuffnorm.tsv"
    message: "cuffnorm is being run"
    params:
        outdir='anl/cuffnorm',
        method='classic-fpkm',
        outformat='simple-table',
        labels=",".join(config["sample_names"])
    threads: 8
    conda: "../envs/cufflinks.yml"
    shell:
        "cuffnorm -p {threads}"
        " -o {params.outdir}"
        " --library-norm-method {params.method}"
        " --output-format {params.outformat}"
        " --labels {params.labels}"
        " {config[genemodel_gtf]}"
        " {input}"
        " > {log} 2>&1"

rule parse_cuffnorm:
    input:
        "anl/cuffnorm/genes.fpkm_table",
        "anl/cuffnorm/isoforms.fpkm_table"
    output:
        "anl/table/abundance.cuffnorm.genes.fpkm.all.xls",
        "anl/table/abundance.cuffnorm.genes.fpkm.exp.xls",
        "anl/table/abundance.cuffnorm.genes.stats.exp.xls",
        "anl/table/abundance.cuffnorm.genes.stats.coexp.xls",
        "anl/table/abundance.cuffnorm.isoforms.fpkm.all.xls",
        "anl/table/abundance.cuffnorm.isoforms.fpkm.exp.xls",
        "anl/table/abundance.cuffnorm.isoforms.stats.exp.xls",
        "anl/table/abundance.cuffnorm.isoforms.stats.coexp.xls"
    log: "logs/parse_cuffnorm.log"
    benchmark: "benchmarks/parse_cuffnorm.tsv"
    message: "parse_cuffnorm is being run"
    shell:
        "python workflow/scripts/parser.py cuffnorm"
        " --outprefix anl/table/abundance.cuffnorm"
        " --gene2desc {config[anno_gene2desc]}"
        " --tran2desc {config[anno_tran2desc]}"
        " > {log} 2>&1"

rule cuffdiff:
    input:
        expand("anl/cuffquant/{sample}/abundances.cxb", sample=config["sample_names"])
    output:
        "anl/cuffdiff/{compare}/gene_exp.diff",
        "anl/cuffdiff/{compare}/isoform_exp.diff"
    log: "logs/cuffdiff/{compare}.log"
    benchmark: "benchmarks/cuffdiff/{compare}.tsv"
    message: "cuffdiff is in progress on {wildcards.compare}"
    params:
        control_cxb=lambda wildcards: find_cxb(wildcards.compare, "control"),
        treatment_cxb=lambda wildcards: find_cxb(wildcards.compare, "treatment"),
        labels=lambda wildcards: make_cuffdiff_label(wildcards.compare),
        outdir='anl/cuffdiff/{compare}',
        dispersion_method='pooled',
        norm_method='geometric'
    threads: 8
    conda: "../envs/cufflinks.yml"
    shell:
        "cuffdiff -p {threads}"
        " -o {params.outdir}"
        " --dispersion-method {params.dispersion_method}"
        " --library-norm-method {params.norm_method}"
        " --labels {params.labels}"
        " {config[genemodel_gtf]}"
        " {params.control_cxb}"
        " {params.treatment_cxb}"
        " > {log} 2>&1"

        
rule parse_cuffdiff:
    input:
        expand(["anl/cuffdiff/{compare}/gene_exp.diff", 
                "anl/cuffdiff/{compare}/isoform_exp.diff"],
               compare=config["compare_names"])
    output:
        "anl/table/deg.genes.xls",
        "anl/table/deg.genes.stats.xls",
        "anl/table/deg.isoforms.xls",
        "anl/table/deg.isoforms.stats.xls",
    log: "logs/parse_cuffdiff.log"
    benchmark: "benchmarks/parse_cuffdiff.tsv"
    message: "parse_cuffdiff is being run"
    params:
        log2fc_criteria=config["log2fc_criteria"],
        p_criteria=config["p_criteria"]
    shell:
        "python workflow/scripts/parser.py cuffdiff"
        " --outprefix anl/table/deg"
        " --gene2desc {config[anno_gene2desc]}"
        " --tran2desc {config[anno_tran2desc]}"
        " --log2fc-criteria {params.log2fc_criteria}"
        " --p-criteria {params.p_criteria}"
        " > {log} 2>&1"

