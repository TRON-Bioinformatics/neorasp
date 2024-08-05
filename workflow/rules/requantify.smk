rule generate_context_fa:
    """
    Convert input table to FASTA file
    """
    input:
        context_file = ""
    output:
        context_fa = "results/idx/context.fa"
    threads: 1
    conda:
        'envs/easyquant.yaml'
    shell:
        'bp_quant '
        'csv2fasta '
        '--input_csv {input.context_file} '
        '--output_fasta {output.context_fa} 2>&1 | tee {log} '


rule bowtie_index:
    input:
        context_fa = 
    params:
        prefix = lambda wildcards, output: os.path.dirname(output.bowtie_index[0]) + "/bowtie"
    threads: 4
    conda:
        "../envs/bowtie2.yaml"
    output:
        bowtie_index = multiext(
            "results/idx/bowtie",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2"
        )
    log: "results/logs/bowtie_index.log"
    shell:
        'bowtie2-build '
        '--threads {threads} '
        '{input.context_fa} '
        '{params.prefix} 2>&1 | tee {log}'
        

rule bowtie_align:
    input:
        ubam = "",
        bowtie_index = rules.bowtie_index.output.bowtie_index
    log: "results/logs/{sample}_bowtie.txt"
    params:
        index_prefix = lambda wildcards, input: input.bowtie_index[0].rstrip(".1.bt2")
    output:
        sam = pipe("results/{sample}/alignment/bowtie_Aligned.out.sam.FIFO")
    threads: 4
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x {params.index_prefix} "
        "-b {input.ubam} "
        "-a --end-to-end "
        "--align-paired-reads "
        "-S {output.sam} "
        "2>&1 | tee {log}"

rule remove_chimeric_alignments:
    """
    Bowtie produces chimeric alignments. Remove and continue with counting
    """
    input:
        sam_pipe = "results/{sample}/alignment/bowtie_Aligned.out.sam.FIFO"
    output:
        sam_filtered_pipe = "results/{sample}/alignment/bowtie_Aligned.out.sam.filtered.FIFO"
    shell:
        'sambamba view '
        '-F not chimeric '
        ''

rule requantify:
    input:
        sam = "results/{sample}/alignment/bowtie_Aligned.out.sam",
        context_seq = get_context_tab
    params:
        requant_dir = lambda wildcards, output: os.path.dirname(output.quant),
        distance = config["requantify"].get('distance', 10),
        interval = "--interval_mode" if config["requantify"].get("interval_mode", True) else "",
        mismatches = "--allow_mismatches" if config["requantify"].get("allow_mismatches", False) else ""
    output:
        quant = "results/{sample}/quantification.tsv",
        read_info = "results/{sample}/read_info.tsv.gz"
    log: "results/logs/{sample}_requantify.log"
    threads: 1
    conda:
        "../envs/easyquant.yaml"
    shell:
        "bp_quant count "
        "-i {input.sam} "
        "-t {input.context_seq} "
        "-d {params.distance} "
        "-o {params.requant_dir} "
        "{params.interval} "
        "{params.mismatches} "
        "2>&1 | tee {log}"
