rule fraser:
    """FRASER
    
    Run FRASER on each sample to calculate PSI and Intron Jaccard Index of splice junctions

    input:
        bam (str): Alignment (BAM sorted by coordinated).
        bai (str): Alignment index.

    output:
        psi_table (str): Splice junction annotated with intron usage

    params:
        working_dir (str): Working directory for fraser functions
        min_read (int): Minimum number of reads to count splice junction. Defaults to 5
        mapq_filter (int): Mapping quality filter to discard low quality alignments in calculations. 
            Defaults to 255
    """

    input:
        bam = rules.samtools.output.bam,
        bai = rules.samtools.output.bai,
    params:
        working_dir = 
            lambda wildcards, output: os.path.dirname(output.psi_table),
        min_read = config['fraser'].get('min_read', 5),
        mapq_filter = config['fraser'].get('mapq_filter', 255)
    log:  "results/{sample}/log/fraser.log"
    output:
        psi_table = "results/{sample}/fraser/junctions_psi.tsv"
    threads: 2
    resources:
        mem_mb = 16000
    container:
        config['container'].get('fraser')
    conda:
        '../envs/fraser.yaml'
    script:
        '../scripts/fraser.R'

rule calculate_junction_cpm:
    """CPM calculation

    Rule to normalize raw splice junction counts obtained by STAR
    using count per million (CPM) metric. This rule normalizes
    the uniquely mapped reads (jCPM_uniquely_mapped), multi mapped
    reads (jCPM_multi_mapped) and total (uniquely and multi) mapped 
    reads (jCPM_total_mapped).

    input:
        star_sj (str): Path to STAR high confidence SJ.out.tab

    output:
        star_sj_cpm (str): Path to normalised SJ counts

    """
    input:
        star_sj = rules.star.output.sj
    output:
        star_sj_cpm = "results/{sample}/fetchdata/parsing/sj_out_tab_cpm.tsv"
    log: "results/{sample}/log/sj_cpm.log"
    threads: 1
    resources:
        mem_mb = 4096
    conda: '../envs/python.yaml'
    container:
        config['container'].get('additional_software')
    script:
	    '../scripts/normalize_star_cpm.py'

rule parse_junctions:
    """Parse junctions

    Rule to parse STAR and FRASER results into standardized splice junction format based on genomic
    coordinates.

    input:
        star_sj (str): Path to STAR high confidence SJ.out.tab
        fraser_psi (str): Path to FRASER metric table
        canonical_junctions (str): Path to canonical junction filter file
        star_cpm (str): Path to CPM normalised splice junction counts.

    output:
        parsed_sj (str): Path to parsed junction table containing 
            STAR and FRASER junctions
        removed_junction (str): Path to expressed canonical junctions

    params:
        read_support (int): Minimum number of unique alignments to report junction. 
            Defaults to 5

    """
    input:
        star_sj = rules.star.output.sj,
        fraser_psi = rules.fraser.output.psi_table,
        star_cpm = rules.calculate_junction_cpm.output.star_sj_cpm,
        canonical_junctions = config['reference']['canonical_juncs']
    output:
        parsed_sj = "results/{sample}/fetchdata/parsing/parsed_sj_star_fraser.tsv",
        removed_junction = "results/{sample}/fetchdata/detected_sj_canonical.tsv"
    log: "results/{sample}/log/sj_parsing.log"
    threads: 1
    resources:
        mem_mb = 8000
    container:
        config['container'].get('splice2neo')
    conda:
        '../envs/R.yaml'
    script:
        '../scripts/parse_junctions.R'


rule filter_mappability:
    """Mapability filter

    Rule to filter junctions located in problematic regions. Problematic regions
    are defined by UCSC and ENCODE lists. These regions are defined as problematic,
    as short-read experiments can not be unambigously mapped to these locations.

    input:
        parsed_sj (str): Path to parsed splice junction table.
        encode_regions (str): Path to ENCODE blacklist.
        ucsc_regions (str): Path to UCSC problematic regions.

    output:
        parsed_sj (str): Path to splice junction table passing the filtering.
        failed_sj (str): Path to splice junction table overlapping problematic regions.
    
    """
    input:
        parsed_sj = rules.parse_junctions.output.parsed_sj,
        encode_regions = config['reference']['encode_mapability'],
        ucsc_regions = config['reference']['ucsc_mapability']
    output:
        parsed_sj = "results/{sample}/fetchdata/mappability/sj_passing_mappability.tsv",
        failed_sj = "results/{sample}/fetchdata/mappability/sj_problematic_mappability.tsv"
    threads: 1
    resources:
        mem_mb = 8000
    container:
        config['container'].get('splice2neo')
    conda:
        '../envs/R.yaml'
    log: "results/{sample}/log/mappability_filter.log"
    script:
        '../scripts/filter_mapability.R'

rule filter_reliable_calls:
    """Filter novel junctions based on expression by user

    input:
        annotated_sj (str):  Path to splice junction table.
    output:
        sj_expression (str): Path to splice junction table with relibale calls.
        sj_low_expression (str): Path to splice junction table failing relibale call parameter.

    """
    input:
        annotated_sj = rules.filter_mappability.output.parsed_sj,
    output:
        sj_expression = temp("results/{sample}/fetchdata/splice2neo/reliable_call/sj_reliable_call.tsv"),
        sj_low_expression =  "results/{sample}/fetchdata/splice2neo/reliable_call/sj_fail_reliable_call.tsv"
    params:
        min_junction_usage = config['reliable_calls'].get('min_junction_usage', 0.01),
        min_junction_cpm = config['reliable_calls'].get('min_junction_cpm', 0.1)
    threads: 1
    resources:
        mem_mb = 8000
    container:
        config['container'].get('splice2neo')
    conda:
        '../envs/R.yaml'
    log: 'results/{sample}/log/filter_reliable_calls.log'
    script:
        '../scripts/filter_reliable_calls.R'

rule add_gene_annotation:
    """Annotation

    Predicted splice junctions are annotated with possible transcript,
    gene and HGNC identifiers. Junctions not overlapping any transcript
    feature from the annotation are removed in this step.

    input:
        parsed_sj (str):  Path to splice junction table.
        transcripts (str): Path to RDS object of reference transcripts.
        gene2hgnc (str): Path to gene to HGNC (gene name) mapping.
        tx2gene (str): Path to transcript to gene mapping.
        rmsk (str): Path to RepeatMasker annotation.
    output:
        annotated_sj (str): Path to splice junction table with feature annotation.
        annotated_sj_problematic (str): Path to table with splice junctions
            not overlapping any feature.

    """
    input:
        parsed_sj = rules.filter_reliable_calls.output.sj_expression,
        transcripts = config['reference']['ref_transcripts'],
        tx2gene = config['reference']['tx2gene'],
        gene2hgnc = config['reference']['gene2symbol'],
        rmsk = config['reference']['rmsk']
    output:
        annotated_sj = "results/{sample}/fetchdata/splice2neo/gene_annot/sj_gene_transcript_overlap.tsv",
        annotated_sj_problematic = "results/{sample}/fetchdata/splice2neo/gene_annot/sj_no_transcript_overlap.tsv"
    threads: 4
    resources:
        mem_mb = 20000
    params:
        extra = "",
    container:
        config['container'].get('splice2neo')
    conda:
        '../envs/R.yaml'
    log:  "results/{sample}/log/add_gene_transcript.log"
    script:
        '../scripts/add_gene_annot.R'

rule filter_gene_hgnc:
    """HGNC filter

    Rule to remove splice junctions from highly polymorphic
    gene loci or genes not located on chromosomes.
    Default genes removed in this step are IG^, TCR^, BCR^ and HLA.

    input:
        parsed_sj (str):  Path to splice junction table.
    output:
        sj_passed_gene (str): Path to table with junctions passing the filter.
        sj_excluded_gene (str): Path to table with junctions falling into exclusion regions.
        sj_gene_exclusion_intention (str): Path to table with detailed information
            of calssification.
    params:
        working_dir (str): Dirname of output files.
        organism (str): Name of organism to apply appropriate gene name filter. Can be human or mouse

    """
    input:
        parsed_sj = rules.add_gene_annotation.output.annotated_sj,
    params:
        working_dir = lambda wildcards, output: os.path.dirname(output.sj_passed_gene),
        organism = config.get('organism', 'human'),
    output:
        sj_excluded_gene = "results/{sample}/fetchdata/splice2neo/gene_name_filter/sj_problematic_gene.tsv",
        sj_passed_gene = "results/{sample}/fetchdata/splice2neo/gene_name_filter/sj_pass_gene.tsv",
        sj_gene_exclusion_intention = "results/{sample}/fetchdata/splice2neo/gene_name_filter/gene_exclusion_intention.tsv"
    threads: 1
    resources:
        mem_mb = 8000
    container:
        config['container'].get('additional_software')
    conda:
        '../envs/python.yaml'
    log:  "results/{sample}/log/gene_filtering.log"
    script:
        '../scripts/filter_gene_name.py'


rule add_context_sequence:
    """Add transcript sequence

    Rule to annotate splice junction candidates with possible transcript sequences.

    input:
        parsed_sj (str):  Path to splice junction table.
        transcripts (str): Path to RDS object of reference transcripts.
        genome (str): Path to 2Bit object of reference genome.
    output:
        annotated_sj (str): Path to splice junction table with context sequences.
        annotated_sj_problematic (str): Path to table with removed junctions.

    """
    input:
        parsed_sj = rules.filter_gene_hgnc.output.sj_passed_gene,
        transcripts = config['reference']['ref_transcripts'],
        genome = config['reference']['2bit']
    output:
        annotated_sj = "results/{sample}/fetchdata/splice2neo/cts/sj_annotated_cts.tsv",
    threads: 4
    resources:
        mem_mb = 20000
    params:
        extra = "",
        cts_size = config['requantify'].get('cts_size', 1000)
    container:
        config['container'].get('splice2neo')
    conda:
        '../envs/R.yaml'
    log:  "results/{sample}/log/add_cts.log"
    script:
        '../scripts/add_tx.R'


rule add_transcript_expression:
    """Add transcript/gene expression

    Rule to annotate splice junction candidate with transcript and 
    gene expression determined with salmon on reference transcripts and
    merging TPM values of StringTie transfrags.

    input:
        annotated_sj (str):  Path to splice junction table.
        transcript_expression (str): Path to salmon transcript level quantification.
        gene_expression (str): Path to salmon gene level quantification.
        transfrags_tpm (str): Path to StringTie transfrags TPM table.
    output:
        sj_expression (str): Path to splice junction table with expression estimates added.

    """
    input:
        annotated_sj = rules.add_context_sequence.output.annotated_sj,
        transcript_expression =  'results/{sample}/salmon_bam/quant.sf',
        gene_expression = 'results/{sample}/salmon_bam/quant.genes.sf',
        transfrags_expression = 'results/{sample}/stringtie/junc_to_tpm.tsv'
    output:
        sj_expression = temp("results/{sample}/fetchdata/splice2neo/sj_annotated_expression.tsv")
    threads: 1
    resources:
        mem_mb = 8000
    container:
        config['container'].get('splice2neo')
    conda:
        '../envs/R.yaml'
    log:  "results/{sample}/log/add_expression_estimates.log"
    script:
        '../scripts/add_tpm.R'

