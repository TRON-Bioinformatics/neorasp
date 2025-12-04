rule stringtie:
    """StringTie reference guided assembly

    Rule to run StringTie reference guided transcript assembly
    on STAR aligned BAM files.

    input:
        bam (string): STAR alignment
        gtf (string): Annotation file (GTF)
    output (string): Path to assembled transcripts
    threads (int): Number of CPUs
    resources:
        mem_mb (int): RAM limit

    """
    input: 
        bam = rules.samtools.output.bam,
        gtf = config['reference']['annotation']
    params:
        junction_count = config['stringtie'].get('min_junc_count', 1),
        anchor = config['stringtie'].get('min_junc_anchor', 10)
    threads: 4
    resources:
        mem_mb = 8000
    container:
        config['container'].get('stringtie')
    output:
        transfrags =
            '<results>/stringtie/transfrags.gtf',
    log:
        "<logs>/stringtie.log"
    benchmark:
        '<benchmarks>/stringtie_bench.txt'
    shell:
        'stringtie '
        '{input.bam} '
        '-o {output.transfrags} '
        '-p {threads} '
        '-a {params.anchor} '
        '-j {params.junction_count} '
        '-G {input.gtf} '
        '-l tx &> {log}'

rule gffcompare:
    input:
        transfrags = rules.stringtie.output.transfrags,
        reference = config['reference']['annotation'],
        reference_genome = config['reference']['genome']
    output:
        gff_stats = '<results>/stringtie/{sample}.stats',
        gff_tmap = '<results>/stringtie/{sample}.transfrags.gtf.tmap',
        gff_refmap = '<results>/stringtie/{sample}.transfrags.gtf.refmap',
        gff_annotated = '<results>/stringtie/{sample}.annotated.gtf',
        gff_loci = '<results>/stringtie/{sample}.loci',
        gff_tracking = '<results>/stringtie/{sample}.tracking',
        junc_to_tx = '<results>/stringtie/junc_to_tx.tab'
    params:
        prefix = lambda wildcards, output: os.path.splitext(output.gff_stats)[0]
    threads: 1
    resources:
        mem_mb = 4000
    container:
        config['container'].get('gffcompare')
    log:
        "<logs>/gffcompare.log"
    benchmark:
        '<benchmarks>/gffcompare_bench.txt'
    shell:
        'gffcompare '
        '-r {input.reference} '
        '-s {input.reference_genome} '
        '-j {output.junc_to_tx} '
        '-o {params.prefix} '
        '{input.transfrags} &> {log}'

rule extract_tpm_from_stringtie:
    input:
        gff_tmap = rules.gffcompare.output.gff_tmap,
    output:
        tpm = '<results>/stringtie/transfrags.tpm.tsv'
    container:
        config['container'].get('additional_software')
    threads: 1
    resources:
        mem_mb = 4000
    log:
        "<logs>/extract_tpm_from_stringtie.log"
    benchmark:
        '<benchmarks>/extract_tpm_from_stringtie_bench.txt'
    script:
        '../scripts/stringtie.py'


rule junc_to_tpm:
    input:
        junc_to_tx = rules.gffcompare.output.junc_to_tx,
        tpm = rules.extract_tpm_from_stringtie.output.tpm
    output:
        junc_to_tpm = '<results>/stringtie/junc_to_tpm.tsv'
    container:
        config['container'].get('additional_software')
    threads: 1
    resources:
        mem_mb = 4000
    log:
        "<logs>/junc_to_stringtie_tpm.log"
    benchmark:
        '<benchmarks>/junc_to_tpm_bench.txt'
    script:
        '../scripts/stringtie.py'
