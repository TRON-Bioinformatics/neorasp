rule prepare_requant:
    """
    Prepare splice junction context sequences
    for re-quantification with easyquant. This
    creates the required input table with unique
    cts identifiers and breakpoints.

    input:
        sj (str): Path to splice junctions with context sequences

    output:
        easyquant_table (str): Path to table with context sequences.
        genes_of_interest (str): A file with gene ids.
    """
    input:
        sj = rules.filter_reliable_calls.output.sj_expression
    output:
        easyquant_table = "results/{sample}/easyquant/context_seq.tsv",
        genes_of_interest = "results/{sample}/easyquant/genes_of_interest.txt"
    threads: 1
    resources:
        mem_mb = 8000
    container:
        'docker://tronbioinformatics/splice2neo:0.6.13'
    conda:
        '../envs/R.yaml'
    log:  "results/{sample}/log/prepare_requantification.log"
    script:
        '../scripts/prepare_quant.R'

rule generate_context_fa:
    """Easyquant FASTA

    Convert easyquant input table to FASTA file

    input:
        easyquant_table (str): Path to table with context sequences
    output:
        context_fa (str): Path to FASTA file of context sequences
    """
    input:
        easyquant_table = rules.prepare_requant.output.easyquant_table
    output:
        context_fa = "results/{sample}/easyquant/context.fa"
    threads: 1
    resources:
        mem_mb = 4000
    container:
        'docker://tronbioinformatics/easyquant:0.6.0'
    conda:
        '../envs/easyquant.yaml'
    log:  "results/{sample}/log/bpquant_csv2fasta.log"
    shell:
        'bp_quant '
        'csv2fasta '
        '--input_csv {input.easyquant_table} '
        '--output_fasta {output.context_fa} 2>&1 | tee {log} '


rule bowtie_index:
    """Bowtie2 index

    Build bowtie2 index of predicted transcript context sequences.

    input:
        context_fa (str): Path to FASTA file with context sequences
    output:
        bowtie_index (List[str]): Paths of bowtie2 index files
    params:
        prefix (str): Path where bowtie2 shall be created
    """
    input:
        context_fa = rules.generate_context_fa.output.context_fa
    params:
        prefix = lambda wildcards, output: output.bowtie_index[0].removesuffix('.1.bt2')
    threads: 4
    resources:
        mem_mb = 8000
    container:
        'docker://tronbioinformatics/easyquant:0.6.0'
    conda:
        '../envs/easyquant.yaml'
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
    """Bowtie2 align

    Targeted mapping of RNA-seq reads against context sequences.
    Alignment generated in this step is used in re-quantification
    counting. Reads mates are aligned against the predicted transcript
    sequence using stringent parameters taken from RSEM. Parameters
    ensure that only properly mapped mates are reported.

    `--no-mixed --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.01`

    input:
        r1 (str): Path to forward read.
        r2 (str): Path to reverse read.
        bowtie_index (List[str]): A list with paths to bowtie2 index files. 
    
    output:
        sam (str): Path to unsorted re-quantification alignment in SAM format.
    
    params:
        index_prefix (str): Prefix of bowtie2 index files.
        report_threshold (str): Number of alignments bowtie2 should report. Defaults to 200.
    
    """
    input:
        r1 = "results/{sample}/fastp/{sample}_R1.fastq.gz",
        r2 = "results/{sample}/fastp/{sample}_R2.fastq.gz",
        bowtie_index = rules.bowtie_index.output.bowtie_index
    log: "results/logs/{sample}_bowtie.txt"
    params:
        index_prefix = lambda wildcards, input: input.bowtie_index[0].removesuffix(".1.bt2"),
        report_threshold = '-a' if config["requantify"].get('bowtie_k_threshold', 200) == 'all' \
            else f'-k {config["requantify"].get("bowtie_k_threshold", 200)}',
    output:
        sam = temp("results/{sample}/easyquant/alignment/bowtie_Aligned.out.sam")
    threads: 4
    resources:
        mem_mb = 16000
    container:
        'docker://tronbioinformatics/easyquant:0.6.0'
    conda:
        '../envs/easyquant.yaml'
    log:  "results/{sample}/log/bowtie_align.log"
    shell:
        "bowtie2 "
        "-p {threads} "
        "-x {params.index_prefix} "
        "{params.report_threshold} "
        "--end-to-end "
        "--no-discordant "
        "--no-mixed "
        "--dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.01 "
        "-1 {input.r1} -2 {input.r2} "
        "> {output.sam} "
        "2> {log}"

rule requantify:
    """Re-quantification

    Counting of junction and spanning read-pairs for predicted
    context sequences. Only properly mapped read pairs are
    counted for each transcript sequence in this step. By
    default, only reads that overlap the junction point
    +/- 10bp mismatch free are counted as junction reads.

    input:
        sam (str):  Path to unsorted re-quantification alignment in SAM format.
        context_seq (str): Path to table with context sequences.

    output:
        quant (str): Path to easyquant quantification table.
        read_info (str): Path to easyquant read into table.

    params:
        requant_dir (str): Path to working directory of easyquant.
        distance (int): Number of mismatch-free bp around junction of interest. Defaults to 10
        interval (str): Execution of easyquant in interval mode. Defaults to True.
        mismatches (str): Allow mismatches in junction region. Defaults to False.
    """
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
    container:
        'docker://tronbioinformatics/easyquant:0.6.0'
    conda:
        '../envs/easyquant.yaml'
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
    """Add quant counts

    Merge re-quantification results back to splice junction table.

    input:
        sj (str):  Path to splice junctions with context sequences
        quantification (str): Path to easyquant quantification table.
    output:
        requantified_sj (str): Path to splice junctions table with merged re-quantification results.
    params:
        requant_dir (str):  Path to working directory of easyquant.
    """
    input:
        sj = rules.filter_reliable_calls.output.sj_expression,
        quantification = rules.requantify.output.quant
    output:
        requantified_sj = temp("results/{sample}/fetchdata/splice2neo/sj_annotated_requantified.tsv")
    params:
        requant_dir = lambda wildcards, input: os.path.dirname(input.quantification)
    log:  "results/{sample}/log/add_requantification_counts.log"
    threads: 1
    resources:
        mem_mb = 8000
    container:
        'docker://tronbioinformatics/splice2neo:0.6.13'
    conda:
        '../envs/R.yaml'
    script:
        '../scripts/quant.R'

rule translate_to_peptide:
    """Peptide sequences

    Annotate splice junctions with mutated peptide sequence.
    In this step, all junctions that generate a mutated coding
    sequence are formated for annotatin with NeoFox.

    input:
        sj (str):  Path to splice junction table.
        cds (str): Path to RDS object of reference coding sequence.
        genome (str): Path to 2Bit object of reference genome.
    output:
        junctions (str): Path to splice junctions table with peptide annotation.
        neofox_annotation (str): Splice junction derived peptides in NeoFox format.

    """
    input:
        sj = rules.add_quant_counts.output.requantified_sj,
        cds = config['reference']['ref_cds'],
        genome = config['reference']['2bit']
    output:
        junctions = "results/{sample}/fetchdata/sj_final.tsv",
        neofox_annotation = "results/{sample}/fetchdata/sj_final_neofox_annotation.tsv"
    log:  "results/{sample}/log/add_peptide_annotation.log"
    threads: 4
    resources:
        mem_mb = 20000
    container:
        'docker://tronbioinformatics/splice2neo:0.6.13'
    conda:
        '../envs/R.yaml'
    script:
        '../scripts/peptide.R'