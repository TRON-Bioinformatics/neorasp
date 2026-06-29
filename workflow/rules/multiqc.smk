rule multiqc:
    """MultiQC
    
    Run multiqc on results to generate a summary report.

    Input:
        get_multiqc_input (List[str]): List of input files for multiqc to summarize.
    Output:
        "results/report/multiqc.html" (str): Path to the generated multiqc report
        directory("results/report/multiqc_data") (str): Directory containing multiqc data files
    params:
        extra (str): Additional command-line arguments for multiqc
        out_dir (str): Output directory for multiqc results
        file_name (str): Name of the multiqc report file
"""
    input:
        get_multiqc_input,
    output:
        "results/report/multiqc.html",
        directory("results/report/multiqc_data"),
    log:
        "results/report/multiqc.log",
    benchmark:
        "results/report/multiqc_bench.txt"
    container:
        config["container"].get("multiqc")
    threads: 1
    params:
        extra="--verbose",
        out_dir=lambda wildcards, output: Path(output[0]).parent,
        file_name=lambda wildcards, output: Path(output[0]).with_suffix("").name,
    shell:
        "multiqc "
        "{params.extra} "
        "--outdir {params.out_dir} "
        "--filename {params.file_name} "
        "{input} "
        "&> {log}"
