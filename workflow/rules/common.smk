#!/usr/bin/env python

def get_final_output():
    """
    Gather final results
    """
    final_files = []
    for sample in samples.itertuples():
        final_files.extend(
            expand(
                "results/{sample}/star/Aligned.sortedByCoord.out.bam", sample=sample.sample_name
            )
        )
        final_files.extend(
            expand(
                "results/{sample}/salmon_bam/quant.sf", sample=sample.sample_name
            )
        )
        final_files.extend(
            expand(
                "results/{sample}/fetchdata/annotated_sj_expression.tsv", sample=sample.sample_name
            )
        )
    return final_files

#rule add_junction_tag:
#    """Add ZJ tag to bam
#     
#    """
#    input:
#        bam = rules.tronmake_expression_star.output.alignment
#    threads: 1
#    resources:
#        mem_mb = 8000
#    params:
#        exe = workflow.source_path('../tools/annotate_ubam_reads.py')
#    log:
#        log = 'results/predictions/star/{sample}/{sample}_bam_annotation.log'
#    conda:
#        '../envs/base.yaml'
#    output:
#        annotated_bam = 'results/predictions/star/{sample}/{sample}_Aligned.annotated.out.bam
#    shell:
#        'python {params.exe} '
#        '--bam {input.bam} '
#        '--cigar True '
#        '> {output.annotate_bam} '
#
#
#rule filter_alignment:
#    """
#    Filter for reads in protein coding regions and unmapped reads
#    """
#
#
#rule revert_bam:
#    """Revert and name sort bam file
#    
#    Rule to revert alignment from STAR into unmapped BAM (uBAM).
#    Attributes are cleared for requantification mapping with
#    easyquant. STAR's special junction tags showing 0-based
#    coordinates of intron and read group string are preserved.
#
#    input:
#        bam (string): Alignment (BAM sorted by name or coodrinate)
#    output:
#        unmapped_bam (string): RNA-seq reads as unmapped bam file
#     
#    """
#    input:
#        bam = rules.tronmake_expression_star.output.alignment,
#    params:
#        tmp_dir = lambda wildcards, output: os.path.join(
#            os.path.dirname(output.unmapped_bam), 'RevertSam_tmp')
#    threads: 1
#    resources:
#        mem_mb = 8000
#    message: 'Creating unmapped bam from STAR bam file...'
#    log:
#        'results/log/star/{sample}_picard_revertSam.log'
#    conda:
#        '../envs/picard.yaml'
#    output:
#        unmapped_bam = 'results/predictions/star/{sample}/{sample}_unmapped.bam'
#    shell:
#        'picard RevertSam '
#        '--INPUT {input.bam} '
#        '--OUTPUT {output.unmapped_bam} '
#        '--ATTRIBUTE_TO_CLEAR NH '
#        '--ATTRIBUTE_TO_CLEAR HI '
#        '--ATTRIBUTE_TO_CLEAR NM '
#        '--ATTRIBUTE_TO_CLEAR MD '
#        '--ATTRIBUTE_TO_CLEAR AS '
#        '--ATTRIBUTE_TO_CLEAR nM '
#        '--ATTRIBUTE_TO_CLEAR jM '
#        '--SANITIZE true '
#        '--KEEP_FIRST_DUPLICATE true '
#        '--TMP_DIR {params.tmp_dir} '
#        '--SORT_ORDER queryname 2>&1 | tee {log} '
#