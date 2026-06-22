rule find_cts_matching_wt:
    input:
        requantified_sj="results/{sample}/fetchdata/splice2neo/sj_annotated_requantified.tsv",
        transcripts=os.path.join(config["index_dir"], "ref_cdna.fa"),
    output:
        blasted_junctions="results/{sample}/fetchdata/blast/requantified_blast_sj.tsv",
    log:
        "results/{sample}/log/blast.log",
    conda:
        "../envs/blast.yaml"
    threads: 8
    resources:
        mem_mb=8000,
    params:
        tmp_dir=lambda wildcards, output: os.path.join(
            os.path.dirname(output[0]), "blast_tmp"
        ),
        exe=workflow.source_path("../scripts/identify_fp_requantification.py"),
    shell:
        "mkdir -p {params.tmp_dir}; "
        "python {params.exe} "
        "--input_splice {input.requantified_sj} "
        "--output {output} "
        "--transcriptome {input.transcripts} "
        "--temp_dir {params.tmp_dir} "
        "--threads {threads} 2>&1 | tee {log}"
