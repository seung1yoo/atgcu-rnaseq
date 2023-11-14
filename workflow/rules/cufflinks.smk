
rule cufflinks:
    input:
        sortBam="anl/star/{sample}/aligned.bam",
    output:
        fpkmGene="anl/cufflinks/{sample}/genes.fpkm_tracking",
        fpkmTran="anl/cufflinks/{sample}/isoforms.fpkm_tracking"
    params:
        outdir='anl/cufflinks/{sample}',
        libraryType=lambda wildcards: config["samples"][wildcards.sample]["cuff_librarytype"]
    threads: 8
    message: "Running Cufflinks on {wildcards.sample}"
    benchmark: "benchmarks/run_cufflinks.{sample}.txt"
    log: 
        stdout="logs/cufflinks/{sample}.stdout",
        stderr="logs/cufflinks/{sample}.stderr",
    conda: "envs/cufflinks.yml"
    shell:
        "cufflinks -p {threads}"
        " --multi-read-correct"
        " --frag-bias-correct {config[genome_fasta]}"
        " --library-type {params.libraryType}"
        " -G {config[genemodel_gtf]}"
        " -o {params.outdir}"
        " {input.sortBam}"

#rule merge_cufflinks:
#    input:
#        fpkm_gene=expand("analysis/cufflinks/{sample_name}/genes.fpkm_tracking", sample_name=config['ordered_sample_list']),
#        fpkm_tran=expand("analysis/cufflinks/{sample_name}/isoforms.fpkm_tracking", sample_name=config['ordered_sample_list'])
#    output:
#        tsv_gene="analysis/cufflinks/genes.fpkm.tsv",
#        tsv_tran="analysis/cufflinks/transcripts.fpkm.tsv"
#    wildcard_constraints:
#        sample_name="[^/]+"
#    params:
#        sample_names=expand("{sample_name}", sample_name=config['ordered_sample_list'])
#    message: "Merge cufflinks outputs"
#    benchmark:
#        "benchmarks/merge_cufflinks.txt"
#    shell:
#        "python {config[pipedir]}/modules_quant/scripts/parser.py"
#        " merge_cufflinks"
#        " {config[workdir]}"
#        " --sample-id-s {params.sample_names}"
#

rule cufflinks:
    input:
        sortBam="anl/star/{sample}/aligned.bam",
    output:
        fpkmGene="anl/cufflinks/{sample}/genes.fpkm_tracking",
        fpkmTran="anl/cufflinks/{sample}/isoforms.fpkm_tracking"
    params:
        outdir='anl/cufflinks/{sample}',
        libraryType=lambda wildcards: config["samples"][wildcards.sample]["cuff_librarytype"]
    threads: 8
    message: "Running Cufflinks on {wildcards.sample}"
    benchmark: "benchmarks/run_cufflinks.{sample}.txt"
    log: 
        stdout="logs/cufflinks/{sample}.stdout",
        stderr="logs/cufflinks/{sample}.stderr",
    conda: "envs/cufflinks.yml"
    shell:
        "cufflinks -p {threads}"
        " --multi-read-correct"
        " --frag-bias-correct {config[genome_fasta]}"
        " --library-type {params.libraryType}"
        " -G {config[genemodel_gtf]}"
        " -o {params.outdir}"
        " {input.sortBam}"
