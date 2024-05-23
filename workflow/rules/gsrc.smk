
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

rule get_mapg2t_gtf:
    input:
        gtf=config["genemodel_gtf"]
    output:
        config["map_gene2tran"]
    log: "logs/get_mapg2t_gtf.log"
    params:
        wkdir="src"
    shell:
        "python workflow/scripts/parser.py" 
        " --wkdir {params.wkdir}"
        " get_mapg2t_gtf"
        " --gtf {input.gtf}"
        " > {log} 2>&1"

rule get_desc:
    input:
        gtf=config["genemodel_gtf"],
        mysqlgene=config["ens_mysql_gene"],
        mysqltran=config["ens_mysql_tran"],
    output:
        config["anno_gene2desc"],
        config["anno_tran2desc"]
    log: "logs/get_desc.log"
    params:
        wkdir="src"
    shell:
        "python workflow/scripts/parser.py" 
        " --wkdir {params.wkdir}"
        " get_desc"
        " --gtf {input.gtf}"
        " --mysqlgene {input.mysqlgene}"
        " --mysqltran {input.mysqltran}"
        " > {log} 2>&1"

