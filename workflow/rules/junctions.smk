rule fraser:
    input:
        bam = rules.tronmake_expression_star.output.alignment,
        gtf = os.path.join(config['index_dir'], 'ref_annot.gtf'),
        strand = rules.tronmake_expression_determine_strandedness.output.strand
    params:
        working_dir = 
            lambda wildcards, output: os.path.dirname(output.psi_table),
        min_read = config['fraser'].get('min_read', 5),
        mapq_filter = config['fraser'].get('mapq_filter', 255),
        exe = workflow.source_path('../scripts/fraser.R')
    log:  "results/{sample}/log/fraser.log"
    output:
        psi_table = "results/{sample}/fraser/junctions_psi.tsv"
    threads: 4
    conda: '../envs/fraser.yaml'
    shell:
        'source {input.strand}; '
        'Rscript --vanilla {params.exe} '
        '--bam {input.bam} '
        '--threads {threads} '
        '--output_table {output.psi_table} '
        '--strandedness ${{featurecounts}} '
        '--min_expression {params.min_read} '
        '--mapq {params.mapq_filter} 2>&1 | tee {log}'

rule parse_junctions:
    input:
        star_sj = rules.tronmake_expression_star.output.sj,
        fraser_psi = rules.fraser.output.psi_table
    output:
        parsed_sj_tmp = "results/{sample}/fetchdata/parsed_sj.tsv.tmp"
    params:
        exe = workflow.source_path('../scripts/parse_junctions.R'),
        read_support = config['fraser'].get('min_read', 5)
    log: "results/{sample}/log/sj_parsing.log"
    threads: 1
    conda: '../envs/R.yaml'
    shell:
        'Rscript --vanilla {params.exe} '
        '--sj {input.star_sj} '
        '--fraser {input.fraser_psi} '
        '--read_support {params.read_support} '
        '--output {output.parsed_sj_tmp} 2>&1 | tee {log}'

rule filter_mapability:
    input:
        parsed_sj = rules.parse_junctions.output.parsed_sj_tmp,
        encode_regions = os.path.join(config['index_dir'], 'mapability', 'encode_blacklist.bed'),
        ucsc_regions =   os.path.join(config['index_dir'], 'mapability', 'ucsc_unusal.bed')
    output:
        parsed_sj = "results/{sample}/fetchdata/parsed_sj.tsv",
        failed_sj = "results/{sample}/fetchdata/parsed_sj_problematic_mapability.tsv"
    threads: 1
    params:
        exe =  workflow.source_path('../scripts/filter_mapability.R')
    conda: '../envs/R.yaml'
    log:  "results/{sample}/log/mapability_filter.log"
    shell:
        'Rscript --vanilla {params.exe} '
         '--juncs {input.parsed_sj} '
         '--encode_blacklist {input.encode_regions} '
         '--ucsc_unusual {input.ucsc_regions} '
         '--output {output.parsed_sj} '
         '--removed_output {output.failed_sj} 2>&1 | tee {log} '

rule add_context_sequence:
    input:
        parsed_sj = rules.filter_mapability.output.parsed_sj,
        transcripts = os.path.join(config['index_dir'], 'ref_transcripts.RDS'),
        tx2gene = os.path.join(config['index_dir'], 'tx2gene.tsv'),
        gene2hgnc = os.path.join(config['index_dir'], 'hgnc_ensembl_id_gencode46.gz'),
        gene_exclusion = os.path.join(config['index_dir'], 'exclusion_pattern.tsv')
    output:
        annotated_sj = "results/{sample}/fetchdata/annotated_sj.tsv",
        annotated_sj_problematic = "results/{sample}/fetchdata/annotated_sj_problematic_gene.tsv"
    params:
        exe = workflow.source_path('../scripts/add_tx.R')
    conda: '../envs/R.yaml'
    log:  "results/{sample}/log/add_cts.log"
    shell:
        'Rscript --vanilla {params.exe} '
        '--juncs {input.parsed_sj} '
        '--transcripts {input.transcripts} '
        '--tx2gene {input.tx2gene} '
        '--gene_exclusion {input.gene_exclusion} '
        '--gene2hgnc {input.gene2hgnc} '
        '--output {output.annotated_sj} '
        '--removed_output {output.annotated_sj_problematic} 2>&1 | tee {log}'

rule add_transcript_expression:
    input:
        annotated_sj = rules.add_context_sequence.output.annotated_sj,
        transcript_expression =  'results/{sample}/salmon_bam/quant.sf',
        gene_expression = 'results/{sample}/salmon_bam/quant.genes.sf'
    output:
        annotated_sj_expression = "results/{sample}/fetchdata/annotated_sj_expression.tsv"
    params:
        exe = workflow.source_path('../scripts/add_tpm.R')
    threads: 1
    conda: '../envs/R.yaml'
    log:  "results/{sample}/log/add_expression_estimates.log"
    shell:
        'Rscript --vanilla {params.exe} '
        '--sj {input.annotated_sj} '
        '--txp {input.transcript_expression} '
        '--gxp {input.gene_expression} '
        '--output {output.annotated_sj_expression} 2>&1 | tee {log}'

#
#rule ctat_splicing:
#    input:
#        sj_out = "results/{sample}/star/SJ.out.tab",
#        chim_junc = "results/{sample}/star/Chimeric.out.junction"
#    output:
#        introns = "results/{sample}/ctat_splicing/{sample}.introns",
#        cancer_introns = "results/{sample}/ctat_splicing/{sample}.cancer.introns"
#    params:
#        prefix = lambda wildcards, output: output.introns.rstrip(".introns"),
#        genome_lib = config['ctat'].get('genome_lib')
#    shell:
#        'python STAR_to_cancer_introns.py '
#        '--ctat_genome_lib {params.genome_lib} '
#        '--SJ_tab_file {input.sj_out} '
#        '--chimJ_file {input.chim_junc} '
#        '--output_prefix {params.prefix } 2>&1 | tee {log}'

