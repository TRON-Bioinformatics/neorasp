rule multiqc:
    input:
        get_multiqc_input
    output:
        "results/report/multiqc.html",
        directory("results/report/multiqc_data"),
    threads: 1
    container:
        config['container'].get("multiqc")
    log:
        "results/report/multiqc.log"
    params:
        extra="--verbose",
        out_dir = lambda wildcards, output:
            Path(output[0]).parent,
        file_name = lambda wildcards, output: 
            Path(output[0]).with_suffix("").name,
    shell:
        "multiqc "
        "{params.extra} "
        "--outdir {params.out_dir} "
        "--filename {params.file_name} "
        "{input} "
        "&> {log}"