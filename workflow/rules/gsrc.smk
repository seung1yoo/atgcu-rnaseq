
#---------- building genomic resources ----------

rule get_anno_uniprot:
    input:
        gtf=config["genemodel_gtf"],
        tsv=config["ens_tsv_uniprot"],
        dat=config["uniprot_dat"]
    output:
        config["anno_gene2go"],
        config["anno_gene2ipr"],
        config["anno_gene2pfam"],
        config["anno_tran2go"],
        config["anno_tran2ipr"],
        config["anno_tran2pfam"],
    log: "logs/get_anno_uniprot.log"
    params:
        wkdir="src"
    shell:
        "python workflow/scripts/parser.py" 
        " --wkdir {params.wkdir}"
        " get_anno_uniprot"
        " --gtf {input.gtf}"
        " --tsv {input.tsv}"
        " --dat {input.dat}"
        " > {log} 2>&1"

