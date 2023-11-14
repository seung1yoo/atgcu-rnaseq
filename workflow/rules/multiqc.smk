



rule multiqc:
    input:
        expand(["anl/fastp/{sample}/{sample}.fastp.json",
                "anl/star/{sample}/Log.final.out",
                "anl/tophat/{sample}/align_summary.txt",
                "anl/hisat2/{sample}/{sample}.stats"],
                sample=config["sample_names"])
    output:
        "anl/multiqc/multiqc_report.html",
        directory("anl/multiqc/multiqc_data"),
    log: "logs/multiqc.log",
    conda: "envs/multiqc.yml"
    shell:
        "multiqc"
        " --outdir anl/multiqc"
        " {input}"


