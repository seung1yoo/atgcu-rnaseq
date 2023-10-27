


rule multiqc_dir:
    input:
        expand("anl/star/{sample}/Log.final.out", sample=config["sorted_sample_names"]),
    output:
        "anl/multiqc/multiqc.html",
        directory("anl/multiqc/multiqc_data"),
    params:
        extra="--data-dir",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "master/bio/multiqc"

