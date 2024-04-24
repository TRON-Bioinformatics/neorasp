rule fraser:
    input:
        bam = "results/hisat2/{sample}/{sample}.sorted.bam",
        gtf = get_genome_annotation
    params:
        working_dir = 
            lambda wildcards, output: os.path.dirname(output.psi_table),
        min_read = config['fraser'].get('min_read', 5),
        mapq_filter = config['fraser'].get('mapq_filter', 60),
    output:
        psi_table = "results/fraser/{sample}/{sample}_junctions.tsv"
    threads: 8
    script:
        '../scripts/fraser.R'

rule regtools:
    """Rule to run Regtools junctions extract/annotate

    Rule to run Regtools junctions subcommand for each sample.
    Regtools extract the expressed splice junctions from the 
    STAR alignment and annotates them with the provided GTF annotation.
    This rule extracts all splice junctions found in the sample and
    not just the junctions in proximity of a variant.

    input:
        bam (string): Alignment (BAM sorted by coord.)
        bai (string): Index of sorted BAM file
        gtf (string): Reference annotation file (GTF)
        fasta (string): Genome fasta file
    output:
        observedJunctions (string): Regtools splice junction table

    """
    input:
        bam = "results/alignment/{sample}/{sample}.cram",
        gtf = get_genome_annotation,
        fasta =  get_genome_fasta
    params:
        anchor_length = config['regtools'].get('anchor', 8),
        strandness = config['regtools'].get('strandness', 0)
    conda:
        '../envs/regtools.yaml'
    log:
        "logs/regtools/{sample}.log"
    output:
        observed_junctions =
            'results/regtools/{sample}/{sample}_expressed_junctions.annot.tdt'
    threads: 1
    shell:
        'regtools junctions extract -a {params.anchor_length} -s {params.strandness} {input.bam} | '
        'regtools junctions annotate -o {output.observed_junctions} - {input.fasta} {input.gtf} &> {log}'


"""
rule targeted_assembly:
    pass

rule splice2neo:
    pass

rule add_gene_expression:
    pass

rule requant:
    pass
""
