rule prepare_requant:
    """Prepare easyquant input
    """
    input:
        sj = rules.add_transcript_expression.output.sj_expression
    output:
        easyquant_table = "results/{sample}/easyquant/context_seq.tsv",
        genes_of_interest = "results/{sample}/easyquant/genes_of_interest.txt"
    params:
        exe = workflow.source_path('../scripts/prepare_quant.R')
    threads: 1
    resources:
        mem_mb = 8000
    conda: '../envs/R.yaml'
    container:
        'docker://tronbioinformatics/splice2neo:0.6.11'
    log:  "results/{sample}/log/prepare_requantification.log"
    shell:
        'Rscript --vanilla {params.exe} '
        '--sj {input.sj} '
        '--output {output.easyquant_table} '
        '--output_genes {output.genes_of_interest} 2>&1 | tee {log}'

rule generate_context_fa:
    """
    Convert input table to FASTA file
    """
    input:
        easyquant_table = rules.prepare_requant.output.easyquant_table
    output:
        context_fa = "results/{sample}/easyquant/context.fa"
    threads: 1
    resources:
        mem_mb = 4000
    conda: '../envs/easyquant.yaml'
    container:
        'docker://tronbioinformatics/easyquant:0.6.0'
    log:  "results/{sample}/log/bpquant_csv2fasta.log"
    shell:
        'bp_quant '
        'csv2fasta '
        '--input_csv {input.easyquant_table} '
        '--output_fasta {output.context_fa} 2>&1 | tee {log} '


rule bowtie_index:
    input:
        context_fa = rules.generate_context_fa.output.context_fa
    params:
        prefix = lambda wildcards, output: os.path.dirname(output.bowtie_index[0]) + "/bowtie"
    threads: 4
    resources:
        mem_mb = 8000
    conda: '../envs/easyquant.yaml'
    container:
        'docker://tronbioinformatics/easyquant:0.6.0'
    log:  "results/{sample}/log/bowtie_index.log"
    output:
        bowtie_index = multiext(
            "results/{sample}/easyquant/idx/bowtie",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2"
        )
    shell:
        'bowtie2-build '
        '--threads {threads} '
        '{input.context_fa} '
        '{params.prefix} 2>&1 | tee {log}'
        

rule bowtie_align:
    input:
        r1 = "results/{sample}/fastp/{sample}_R1.fastq.gz",
        r2 = "results/{sample}/fastp/{sample}_R2.fastq.gz",
        bowtie_index = rules.bowtie_index.output.bowtie_index
    log: "results/logs/{sample}_bowtie.txt"
    params:
        index_prefix = lambda wildcards, input: input.bowtie_index[0].rstrip(".1.bt2"),
        report_threshold = config["requantify"].get('bowtie_k_threshold', 200),
    output:
        sam = temp("results/{sample}/easyquant/alignment/bowtie_Aligned.out.sam")
    threads: 4
    resources:
        mem_mb = 16000
    conda: '../envs/easyquant.yaml'
    container:
        'docker://tronbioinformatics/easyquant:0.6.0'
    log:  "results/{sample}/log/bowtie_align.log"
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x {params.index_prefix} "
        "-k {params.report_threshold} "
        "--end-to-end "
        "--no-discordant "
        "--no-mixed "
        "--dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.01 "
        "-1 {input.r1} -2 {input.r2} "
        "> {output.sam} "
        "2> {log}"

rule requantify:
    input:
        sam = "results/{sample}/easyquant/alignment/bowtie_Aligned.out.sam",
        context_seq = rules.prepare_requant.output.easyquant_table
    params:
        requant_dir = lambda wildcards, output: os.path.dirname(output.quant),
        distance = config["requantify"].get('distance', 10),
        interval = "--interval_mode" if config["requantify"].get("interval_mode", True) else "",
        mismatches = "--allow_mismatches" if config["requantify"].get("allow_mismatches", False) else ""
    output:
        quant = "results/{sample}/easyquant/quantification.tsv",
        read_info = "results/{sample}/easyquant/read_info.tsv.gz"
    log: "results/logs/{sample}_requantify.log"
    threads: 1
    resources:
        mem_mb = 8000
    conda: '../envs/easyquant.yaml'
    container:
        'docker://tronbioinformatics/easyquant:0.6.0'
    log:  "results/{sample}/log/requantification.log"
    shell:
        "bp_quant count "
        "-i {input.sam} "
        "-t {input.context_seq} "
        "-d {params.distance} "
        "-o {params.requant_dir} "
        "{params.interval} "
        "{params.mismatches} "
        "2>&1 | tee {log}"

rule add_quant_counts:
    input:
        sj = rules.add_transcript_expression.output.sj_expression,
        quantification = rules.requantify.output.quant
    output:
        requantified_sj = "results/{sample}/fetchdata/splice2neo/sj_annotated_requantified.tsv"
    params:
        exe = workflow.source_path('../scripts/quant.R'),
        requant_dir = lambda wildcards, input: os.path.dirname(input.quantification)
    log:  "results/{sample}/log/add_requantification_counts.log"
    threads: 1
    resources:
        mem_mb = 8000
    conda: '../envs/R.yaml'
    container:
        'docker://tronbioinformatics/splice2neo:0.6.11'
    shell:
        'Rscript --vanilla {params.exe} '
        '--sj {input.sj} '
        '--requant {params.requant_dir} '
        '--output {output.requantified_sj} 2>&1 | tee {log} '

rule translate_to_peptide:
    input:
        sj = rules.add_quant_counts.output.requantified_sj,
        cds = os.path.join(config['index_dir'], 'ref_cds.RDS'),
        genome = os.path.join(config['index_dir'], 'ref_genome.2bit')
    output:
        junctions = "results/{sample}/fetchdata/sj_final_results.tsv",
        neofox_annotation = "results/{sample}/fetchdata/sj_results_neofox_annotation.tsv"
    params:
        exe = workflow.source_path('../scripts/peptide.R'),
    log:  "results/{sample}/log/add_peptide_annotation.log"
    threads: 1
    resources:
        mem_mb = 16000
    conda: '../envs/R.yaml'
    container:
        'docker://tronbioinformatics/splice2neo:0.6.11'
    shell:
        'Rscript --vanilla {params.exe} '
        '--sj {input.sj} '
        '--sample {wildcards.sample} '
        '--cds {input.cds} '
        '--genome {input.genome} '
        '--output {output.junctions} '
        '--output_neofox {output.neofox_annotation} 2>&1 | tee {log}'