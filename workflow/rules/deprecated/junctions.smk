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
    exe (str): Path to python script

"""
    input:
        parsed_sj=rules.add_gene_annotation.output.annotated_sj,
    output:
        sj_excluded_gene="results/{sample}/fetchdata/splice2neo/sj_problematic_gene.tsv",
        sj_passed_gene=temp("results/{sample}/fetchdata/splice2neo/sj_pass_gene.tsv"),
        sj_gene_exclusion_intention="results/{sample}/fetchdata/splice2neo/gene_exclusion_intention.tsv",
    log:
        "results/{sample}/log/gene_filtering.log",
    conda:
        "../envs/python.yaml"
    container:
        "docker://tronbioinformatics/tron_data_utils:0.0.1"
    threads: 1
    resources:
        mem_mb=8000,
    params:
        working_dir=lambda wildcards, output: os.path.dirname(output.sj_passed_gene),
        exe=workflow.source_path("../scripts/filter_gene_regex.py"),
    shell:
        "python {params.exe} "
        "-i {input.parsed_sj} "
        "-o {params.working_dir} 2>&1 | tee {log}"


# rule bedgraph_to_bigwig:
#    """BigWig creation
#
#    Convert coverage bedGraph files from STAR to binary
#    BigWig for genome wide coverage in IGV.
#
#    input:
#        bedGraph_forward (str): Coverage of forward strand in bedGraph format.
#        bedGraph_reverse (str): Coverage of reverse strand in bedGraph format.
#        chromsizes (str): Path to TSV file describing chromosome sizes.
#    output:
#        bw_forward (str): Coverage of forward strand in BigWig format.
#        bw_reverse (str): Coverage of reverse strand in BigWig format.
#    params:
#        extra (str): Additional parameters passed to bedGraphToBigWig.
#    """
#    input:
#        bedGraph_forward = rules.star.output.forward_wig,
#        bedGraph_reverse = rules.star.output.reverse_wig,
#        chromsizes = config['reference']['chromsizes']
#    output:
#        bw_forward = "results/{sample}/star/Signal.Unique.str1.bw",
#        bw_reverse = "results/{sample}/star/Signal.Unique.str2.bw"
#    log:
#        "results/{sample}/log/bedgraph2bigwig.log"
#    params:
#        extra = "" # optional params string
#    container:
#        'docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:472--h9b8f530_1'
#    conda:
#        '../envs/ucsc_bedgraph_to_bigwig.yaml'
#    shell:
#        '''
#        bedGraphToBigWig {params.extra} {input.bedGraph_forward} {input.chromsizes} {output.bw_forward} &> {log}
#        bedGraphToBigWig {params.extra} {input.bedGraph_reverse} {input.chromsizes} {output.bw_reverse} &>> {log}
#        '''
