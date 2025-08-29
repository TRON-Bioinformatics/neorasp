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
        'docker://quay.io/biocontainers/stringtie:3.0.1--h00789bb_0',
    output:
        transfrags =
            'results/{sample}/stringtie/transfrags.gtf',
    log:
        "results/{sample}/log/stringtie.log"
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
        gff_stats = 'results/{sample}/stringtie/{sample}.stats',
        gff_tmap = 'results/{sample}/stringtie/{sample}.transfrags.gtf.tmap',
        gff_refmap = 'results/{sample}/stringtie/{sample}.transfrags.gtf.refmap',
        gff_annotated = 'results/{sample}/stringtie/{sample}.annotated.gtf',
        gff_loci = 'results/{sample}/stringtie/{sample}.loci',
        gff_tracking = 'results/{sample}/stringtie/{sample}.tracking',
        junc_to_tx = 'results/{sample}/stringtie/junc_to_tx.tab'
    params:
        prefix = lambda wildcards, output: os.path.splitext(output.gff_stats)[0]
    threads: 1
    container:
        'docker://quay.io/biocontainers/gffcompare:0.12.10--h9948957_0'
    log:
        "results/{sample}/log/gffcompare.log"
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
        tpm = 'results/{sample}/stringtie/transfrags.tpm.tsv'
    run:
        import pandas as pd
        df = pd.read_csv(input.gff_tmap, sep='\t')
        df = df[['qry_id', 'TPM']]
        df.columns = ['stringtie_tx_id', 'stringtie_TPM']
        df.to_csv(output.tpm, sep='\t', index=False)

rule junc_to_tpm:
    input:
        junc_to_tx = rules.gffcompare.output.junc_to_tx,
        tpm = rules.extract_tpm_from_stringtie.output.tpm
    output:
        junc_to_tpm = 'results/{sample}/stringtie/junc_to_tpm.tsv'
    run:
        import pandas as pd
        junc_df = pd.read_csv(input.junc_to_tx, sep='\t', header=None)
        junc_df.columns = ['chrom', 'start', 'end', 'strand', 'tx_ids']
        tpm_df = pd.read_csv(input.tpm, sep='\t')
        junc_df['tx_ids'] = junc_df['tx_ids'].str.split(',')
        junc_df = junc_df.explode('tx_ids').reset_index(drop=True)
        merged_df = pd.merge(junc_df, tpm_df, left_on='tx_ids', right_on='stringtie_tx_id', how='left')
        merged_df = merged_df[['chrom', 'start', 'end', 'strand', 'tx_ids', 'stringtie_TPM']]
        merged_df = merged_df.groupby(['chrom', 'start', 'end', 'strand']).agg({'tx_ids': lambda x: ','.join(x), 'stringtie_TPM': 'sum'}).reset_index()
        merged_df.to_csv(output.junc_to_tpm, sep='\t', index=False)

#rule gtf_to_fasta:
#    input:
#        gtf=rules.stringtie.output.transfrags
#        ref= config
#    output:
#        fasta="results/{sample}/stringtie/transfrags.fa"
#    shell:
#        """
#        gffread {input.gtf} -g {input.ref} -w {output.fasta}
#        """

#rule run_rnasamba:
#    input:
#        fasta=rules.gtf_to_fasta.output.fasta
#		model=config
#    output:
#        protein_tsv = "results/{sample}/stringtie/transfrags.rnasamba.tsv",
#        protein_fasta = "results/{sample}/stringtie/transfrags.rnasamba.protein.fasta"
#	shell:
#        """
#        rnasamba classify -p {output.protein_fasta} {output.protein_tsv} {input.fasta} {input.model}
#        """
