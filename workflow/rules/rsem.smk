
def is_exists_rsem_index():
    prefix = config["rsem_index"]
    for _ext in ["seq", "grp", "ti", ".transcripts.fa", ".idx.fa", ".n2g.idx.fa"]:
        if not Path(f"{prefix}.{_ext}").exists():
            return 1
    return 0

if is_exists_rsem_index():
    rule rsem_index:
        input:
            reference_genome=config["genome_fasta"],
        output:
            seq="src/rsem_index/reference.seq",
            grp="src/rsem_index/reference.grp",
            ti="src/rsem_index/reference.ti",
            tran_fa="src/rsem_index/reference.transcripts.fa",
            idx_fa="src/rsem_index/reference.idx.fa",
            n2g_idx_fa="src/rsem_index/reference.n2g.idx.fa",
        params:
            extra=f"--gtf {config['genemodel_gtf']}",
        log: "logs/rsem_index.log",
        benchmark: "benchmarks/rsem_index.tsv"
        wrapper:
            "v3.10.2/bio/rsem/prepare-reference"


rule rsem_calculate_expression:
    input:
        bam=rules.star_align.output.tranBam,
        #fq_one="anl/fastp/{sample}/{sample}_R1.fastq.gz",
        #fq_two="anl/fastp/{sample}/{sample}_R2.fastq.gz",
        reference=multiext("src/rsem_index/reference",
                ".seq", 
                ".grp", 
                ".ti", 
                ".transcripts.fa", 
                ".idx.fa", 
                ".n2g.idx.fa"
                ),
    output:
        genes_results="anl/rsem/{sample}/rsem.genes.results",
        isoforms_results="anl/rsem/{sample}/rsem.isoforms.results",
    params:
        paired_end=True,
        extra="--seed 42",
    log: "logs/rsem_calculate_expression/{sample}.log",
    benchmark: "benchmarks/rsem_calculate_expression/{sample}.tsv"
    threads: 8
    wrapper:
        "v3.10.2/bio/rsem/calculate-expression"

rule parse_rsem:
    input:
        gene=expand("anl/rsem/{sample}/rsem.genes.results", sample=config['sample_names']),
        tran=expand("anl/rsem/{sample}/rsem.isoforms.results", sample=config['sample_names']),
    output:
        gene_fpkm="anl/table/abundance.rsem.genes.fpkm.all.xls",
        gene_fpkm_exp="anl/table/abundance.rsem.genes.fpkm.exp.xls",
        gene_tpm="anl/table/abundance.rsem.genes.tpm.xls",
        gene_count="anl/table/abundance.rsem.genes.readcount.xls",
        gene_length="anl/table/abundance.rsem.genes.length.xls",
        gene_stat_exp="anl/table/abundance.rsem.genes.stats.exp.xls",
        
        tran_fpkm="anl/table/abundance.rsem.isoforms.fpkm.all.xls",
        tran_fpkm_exp="anl/table/abundance.rsem.isoforms.fpkm.exp.xls",
        tran_tpm="anl/table/abundance.rsem.isoforms.tpm.xls",
        tran_count="anl/table/abundance.rsem.isoforms.readcount.xls",
        tran_length="anl/table/abundance.rsem.isoforms.length.xls",
        tran_stat_exp="anl/table/abundance.rsem.isoforms.stats.exp.xls",

    log: "logs/parse_rsem.log"
    benchmark: "benchmarks/parse_rsem.tsv"
    message: "parse_rsem is being run"
    params:
        samples=expand("{sample}", sample=config['sample_names'])
    shell:
        "python workflow/scripts/parser.py rsem"
        " --outprefix anl/table/abundance.rsem"
        " --gene2desc {config[anno_gene2desc]}"
        " --tran2desc {config[anno_tran2desc]}"
        " --sample-id-s {params.samples}"
        " > {log} 2>&1"

