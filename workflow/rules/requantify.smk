rule prepare_requant:
    """Prepare easyquant input
    """
    input:
        sj = rules.add_transcript_expression.output.sj_expression
    output:
        easyquant_table = "results/{sample}/fetchdata/easyquant/context_seq.tsv",
        genes_of_interest = "results/{sample}/fetchdata/easyquant/genes_of_interest.txt"
    params:
        exe = workflow.source_path('../scripts/prepare_quant.R')
    threads: 1
    resources:
        mem_mb = 8000
    conda: '../envs/R.yaml'
    log:  "results/{sample}/log/prepare_requantification.log"
    shell:
        'Rscript --vanilla {params.exe} '
        '--sj {input.sj} '
        '--output {output.easyquant_table} '
        '--output_genes {output.genes_of_interest} 2>&1 | tee {log}'

rule subset_reads:
    input:
        unmapped_fq1 = rules.tronmake_expression_star.output.unmapped_fq1,
        unmapped_fq2 = rules.tronmake_expression_star.output.unmapped_fq2,
        bam = rules.tronmake_expression_star.output.alignment,
        genes_of_interest = rules.prepare_requant.output.genes_of_interest,
        gtf = os.path.join(config['index_dir'], 'ref_annot.gtf'),
    output:
        unmapped_bam = "results/{sample}/fetchdata/easyquant/target_reads.ubam"
    shadow: "shallow"
    shell:
        '''
        genes_joined=$(cat {input.genes_of_interest})

        grep -v '^##' {input.gtf} | \
            awk '$3 == "gene"' | \
                grep -E "gene_id "($genes_joined)"" | \
                    awk '{{print $1 "\t" $4-1 "\t" $5 "\t" "g" "\t" 0 "\t" $7}}' | bedtools merge -s > gene_regions.bed
        

        samtools view -b -h --fetch-pairs -L gene_regions.bed {input.bam} | \
            samtools collate -u -O | \
                samtools fastq -1 extracted_R1.fastq -2 extracted_R2.fastq -s /dev/null -0 /dev/null

        cat extracted_R1.fastq {input.unmapped_fq1} > target_reads_R1.fastq
        cat extracted_R2.fastq {input.unmapped_fq2} > target_reads_R2.fastq

        samtools import -1 target_reads_R1.fastq -2 target_reads_R2.fastq -o {output.unmapped_bam} 
        '''

rule generate_context_fa:
    """
    Convert input table to FASTA file
    """
    input:
        easyquant_table = rules.prepare_requant.output.easyquant_table
    output:
        context_fa = "results/{sample}/fetchdata/easyquant/context.fa"
    threads: 1
    resources:
        mem_mb = 4000
    conda: '../envs/easyquant.yaml'
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
    log:  "results/{sample}/log/bowtie_index.log"
    output:
        bowtie_index = multiext(
            "results/{sample}/fetchdata/easyquant/idx/bowtie",
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
        unpack(tronmake_expression.get_fq),
        bowtie_index = rules.bowtie_index.output.bowtie_index
    log: "results/logs/{sample}_bowtie.txt"
    params:
        index_prefix = lambda wildcards, input: input.bowtie_index[0].rstrip(".1.bt2"),
        report_threshold = config["requantify"].get('bowtie_k_threshold', 200),
    output:
        sam = temp("results/{sample}/fetchdata/easyquant/alignment/bowtie_Aligned.out.sam")
    threads: 4
    resources:
        mem_mb = 16000
    conda: '../envs/easyquant.yaml'
    log:  "results/{sample}/log/bowtie_align.log"
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x {params.index_prefix} "
        "-k {params.report_theshold} "
        "--end-to-end "
        "--no-discordant "
        "--no-mixed "
        "--no-unal "
        "--dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.02 "
        "-1 {input.r1} -2 {input.r2} "
        "> {output.sam} "
        "2> {log}"

rule requantify:
    input:
        sam = "results/{sample}/fetchdata/easyquant/alignment/bowtie_Aligned.out.sam",
        context_seq = rules.prepare_requant.output.easyquant_table
    params:
        requant_dir = lambda wildcards, output: os.path.dirname(output.quant),
        distance = config["requantify"].get('distance', 10),
        interval = "--interval_mode" if config["requantify"].get("interval_mode", True) else "",
        mismatches = "--allow_mismatches" if config["requantify"].get("allow_mismatches", False) else ""
    output:
        quant = "results/{sample}/fetchdata/easyquant/quantification.tsv",
        read_info = "results/{sample}/fetchdata/easyquant/read_info.tsv.gz"
    log: "results/logs/{sample}_requantify.log"
    threads: 1
    resources:
        mem_mb = 8000
    conda: '../envs/easyquant.yaml'
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
        requantified_sj = "results/{sample}/fetchdata/requantified_sj.tsv"
    params:
        exe = workflow.source_path('../scripts/quant.R'),
        requant_dir = lambda wildcards, input: os.path.dirname(input.quantification)
    log:  "results/{sample}/log/add_requantification_counts.log"
    threads: 1
    resources:
        mem_mb = 8000
    conda: '../envs/R.yaml'
    shell:
        'Rscript --vanilla {params.exe} '
        '--sj {input.sj} '
        '--requant {params.requant_dir} '
        '--output {output.requantified_sj} 2>&1 | tee {log} '
