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
        bam=rules.samtools.output.bam,
        bai=rules.samtools.output.bai,
    output:
        psi_table="<results>/fraser/junctions_psi.tsv",
    log:
        "<logs>/fraser.log",
    benchmark:
        "<benchmarks>/fraser_bench.txt"
    conda:
        "../envs/fraser.yaml"
    container:
        config["container"].get("fraser")
    threads: 2
    resources:
        mem_mb=16000,
    params:
        working_dir=lambda wildcards, output: os.path.dirname(output.psi_table),
        min_read=config["fraser"].get("min_read", 5),
        mapq_filter=config["fraser"].get("mapq_filter", 255),
    script:
        "../scripts/fraser.R"


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
        star_sj=rules.star.output.sj,
    output:
        star_sj_cpm="<results>/fetchdata/parsing/sj_out_tab_cpm.tsv",
    log:
        "<logs>/sj_cpm.log",
    benchmark:
        "<benchmarks>/calculate_junction_cpm_bench.txt"
    conda:
        "../envs/python.yaml"
    container:
        config["container"].get("additional_software")
    threads: 1
    resources:
        mem_mb=4096,
    script:
        "../scripts/normalize_star_cpm.py"


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
        star_sj=rules.star.output.sj,
        fraser_psi=rules.fraser.output.psi_table,
        star_cpm=rules.calculate_junction_cpm.output.star_sj_cpm,
        canonical_junctions=config["reference"]["canonical_juncs"],
    output:
        parsed_sj="<results>/fetchdata/parsing/parsed_sj_star_fraser.tsv",
        removed_junction="<results>/fetchdata/detected_sj_canonical.tsv",
    log:
        "<logs>/sj_parsing.log",
    benchmark:
        "<benchmarks>/parse_junctions_bench.txt"
    conda:
        "../envs/R.yaml"
    container:
        config["container"].get("splice2neo")
    threads: 1
    resources:
        mem_mb=8000,
    script:
        "../scripts/parse_junctions.R"


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
        parsed_sj=rules.parse_junctions.output.parsed_sj,
        encode_regions=config["reference"]["encode_mapability"],
        ucsc_regions=config["reference"]["ucsc_mapability"],
    output:
        parsed_sj="<results>/fetchdata/mappability/sj_passing_mappability.tsv",
        failed_sj="<results>/fetchdata/mappability/sj_problematic_mappability.tsv",
    log:
        "<logs>/mappability_filter.log",
    benchmark:
        "<benchmarks>/filter_mappability_bench.txt"
    conda:
        "../envs/R.yaml"
    container:
        config["container"].get("splice2neo")
    threads: 1
    resources:
        mem_mb=8000,
    script:
        "../scripts/filter_mapability.R"


rule filter_reliable_calls:
    """Filter novel junctions based on expression by user

input:
    annotated_sj (str):  Path to splice junction table.
output:
    sj_expression (str): Path to splice junction table with relibale calls.
    sj_low_expression (str): Path to splice junction table failing relibale call parameter.

"""
    input:
        annotated_sj=rules.filter_mappability.output.parsed_sj,
    output:
        sj_expression=temp(
            "<results>/fetchdata/splice2neo/reliable_call/sj_reliable_call.tsv"
        ),
        sj_low_expression="<results>/fetchdata/splice2neo/reliable_call/sj_fail_reliable_call.tsv",
    log:
        "<logs>/filter_reliable_calls.log",
    benchmark:
        "<benchmarks>/filter_reliable_calls_bench.txt"
    conda:
        "../envs/R.yaml"
    container:
        config["container"].get("splice2neo")
    threads: 1
    resources:
        mem_mb=8000,
    params:
        min_junction_usage=config["reliable_calls"].get("min_junction_usage", 0.01),
        min_junction_cpm=config["reliable_calls"].get("min_junction_cpm", 0.1),
    script:
        "../scripts/filter_reliable_calls.R"


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
        parsed_sj=rules.filter_reliable_calls.output.sj_expression,
        transcripts=config["reference"]["ref_transcripts"],
        cds=config["reference"]["ref_cds"],
        tx2gene=config["reference"]["tx2gene"],
        gene2hgnc=config["reference"]["gene2symbol"],
        rmsk=config["reference"]["rmsk"],
    output:
        annotated_sj="<results>/fetchdata/splice2neo/gene_annot/sj_gene_transcript_overlap.tsv",
        annotated_sj_problematic="<results>/fetchdata/splice2neo/gene_annot/sj_no_transcript_overlap.tsv",
        annotated_sj_non_coding="<results>/fetchdata/splice2neo/gene_annot/sj_nc_gene_transcript_overlap.tsv",
    log:
        "<logs>/add_gene_transcript.log",
    benchmark:
        "<benchmarks>/add_gene_annotation_bench.txt"
    conda:
        "../envs/R.yaml"
    container:
        config["container"].get("splice2neo")
    threads: 4
    resources:
        mem_mb=20000,
    params:
        extra="",
    script:
        "../scripts/add_gene_annot.R"


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
        parsed_sj=rules.add_gene_annotation.output.annotated_sj,
    output:
        sj_excluded_gene="<results>/fetchdata/splice2neo/gene_name_filter/sj_problematic_gene.tsv",
        sj_passed_gene="<results>/fetchdata/splice2neo/gene_name_filter/sj_pass_gene.tsv",
        sj_gene_exclusion_intention="<results>/fetchdata/splice2neo/gene_name_filter/gene_exclusion_intention.tsv",
    log:
        "<logs>/gene_filtering.log",
    benchmark:
        "<benchmarks>/filter_gene_hgnc_bench.txt"
    conda:
        "../envs/python.yaml"
    container:
        config["container"].get("additional_software")
    threads: 1
    resources:
        mem_mb=8000,
    params:
        working_dir=lambda wildcards, output: os.path.dirname(output.sj_passed_gene),
        organism=config.get("organism", "human"),
    script:
        "../scripts/filter_gene_name.py"


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
        annotated_sj="<results>/fetchdata/splice2neo/gene_name_filter/sj_pass_gene.tsv",
        transcript_expression="<results>/salmon_bam/quant.sf",
        gene_expression="<results>/salmon_bam/quant.genes.sf",
        transfrags_expression="<results>/stringtie/junc_to_tpm.tsv",
    output:
        sj_expression=temp("<results>/fetchdata/splice2neo/sj_annotated_expression.tsv"),
    log:
        "<logs>/add_expression_estimates.log",
    benchmark:
        "<benchmarks>/add_transcript_expression_bench.txt"
    conda:
        "../envs/R.yaml"
    container:
        config["container"].get("splice2neo")
    threads: 1
    resources:
        mem_mb=8000,
    script:
        "../scripts/add_tpm.R"


checkpoint split_junc_table:
    """
Split detected splice junctions into chunks of fixed size for
parallel processing with splice2neo. By default, the chunks
will have the prefix `splice2neo_input_` followed by a decimal number.

input:
    parsed_sj (str): Splice junctions table containing all columns required for splice2neo annotation
output:
    A directory for storing chunked TSV files.
params:
    scatter_size (int): Number of junctions to include in chunk. Default is 100.

"""
    input:
        parsed_sj="<results>/fetchdata/splice2neo/sj_annotated_expression.tsv",
    output:
        directory("<results>/fetchdata/splice2neo/split_junc"),
    log:
        "<logs>/split_junc.log",
    benchmark:
        "<benchmarks>/split_junc_table_bench.txt"
    container:
        config["container"].get("shell_utils")
    threads: 1
    resources:
        mem_mb=8000,
    params:
        extra="",
        scatter_size=config["splice2neo"].get("scatter_size", 100),
        prefix=lambda wildcards, output: os.path.join(output[0], "splice2neo_input_"),
    shell:
        """
        mkdir -p {output}
        head -n1 {input.parsed_sj} > {output}/header.tsv
        
        split -d -l {params.scatter_size} <(tail -n+2 {input.parsed_sj}) {params.prefix}
        wait $!
        
        for file in {params.prefix}*
        do
            cat {output}/header.tsv ${{file}} > {output}/temp_file && mv {output}/temp_file ${{file}}
        done
        """


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
        parsed_sj="<results>/fetchdata/splice2neo/split_junc/splice2neo_input_{chunkID}",
        transcripts=config["reference"]["ref_transcripts"],
        #genome = config['reference']['genome']
        genome=config["reference"]["2bit"],
    output:
        annotated_sj="<results>/fetchdata/splice2neo/cts/sj_annotated_cts_{chunkID}.tsv",
    log:
        "<logs>/add_cts_{chunkID}.log",
    benchmark:
        "<benchmarks>/add_context_sequence_{chunkID}_bench.txt"
    conda:
        "../envs/R.yaml"
    container:
        config["container"].get("splice2neo")
    threads: 1
    resources:
        mem_mb=20000,
    params:
        extra="",
        cts_size=config["requantify"].get("cts_size", 1000),
    script:
        "../scripts/add_tx.R"


rule translate_to_peptide:
    """Peptide sequences

Annotate splice junctions with mutated peptide sequence.
In this step, all junctions that generate a mutated coding
sequence are formated for annotation with NeoFox.

input:
    sj (str):  Path to splice junction table.
    cds (str): Path to RDS object of reference coding sequence.
    genome (str): Path to 2Bit object of reference genome.
output:
    peptide_junc (str): Path to splice junctions table with peptide annotation.
    peptide_fasta (str): Splice junction derived protein sequences in FASTA format e.g. for MassSpec
    neofox_annotation (str): plice junction derived peptides in NeoFox format.

"""
    input:
        annotated_sj="<results>/fetchdata/splice2neo/cts/sj_annotated_cts_{chunkID}.tsv",
        cds=config["reference"]["ref_cds"],
        genome=config["reference"]["2bit"],
    output:
        peptide_junc="<results>/fetchdata/splice2neo/pep/sj_annotated_peptide_{chunkID}.tsv",
        peptide_fasta="<results>/fetchdata/splice2neo/pep/sj_annotated_peptide_{chunkID}.fasta",
        neofox_annotation="<results>/fetchdata/splice2neo/pep/sj_neofox_annotation_{chunkID}.tsv",
    log:
        "<logs>/add_peptide_annotation_{chunkID}.log",
    benchmark:
        "<benchmarks>/translate_to_peptide_{chunkID}_bench.txt"
    conda:
        "../envs/R.yaml"
    container:
        config["container"].get("splice2neo")
    threads: 1
    resources:
        mem_mb=20000,
    params:
        peptide_flank_size=config["splice2neo"].get("peptide_flank_size", 13),
    script:
        "../scripts/peptide2.R"


rule gather_splice2neo:
    """
Gather scattered junction anntotation files from splice2neo
and merge into unified tables for further processing. Input
of this rule are splice2neo annotated chunks.

input:
    peptide_junc (List[str]): Paths to peptide annotated splice junctions.
    peptide_fasta (List[str]): Paths to splice junctions derived proteins in FASTA.
    neofox_annotation (List[str]): Paths to neofox annotation files.

output:
    sj_annot_cts_peptide (str): Path to transcript and peptide annotated splice junction table.
    peptide_fasta (str): Path to peptide FASTA file
    neofox_annotation (str): Path to neofox annotation file

"""
    input:
        unpack(aggregate_splice2neo_output),
    output:
        sj_annot_cts_peptide="<results>/fetchdata/splice2neo/sj_annotated_peptide.tsv",
        peptide_fasta="<results>/fetchdata/sj_final_peptides.fasta",
        neofox_annotation="<results>/fetchdata/sj_final_neofox_annotation.tsv",
    log:
        "<logs>/splice2neo_gather.log",
    benchmark:
        "<benchmarks>/gather_splice2neo_bench.txt"
    conda:
        "../envs/R.yaml"
    container:
        config["container"].get("shell_utils")
    threads: 1
    resources:
        mem_mb=8000,
    script:
        "../scripts/gather_splice2neo.sh"
