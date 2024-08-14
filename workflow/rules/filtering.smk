
#rule select_cancer_junctions:
#    input:
#        requantified_sj = "results/{sample}/fetchdata/blasted_sj.tsv",
#        filter_list = os.path.join(config['index_dir'], ''),
#        filter_list_annotation = os.path.join(config['index_dir'], ''),
#    output:
#        tumor_specific_junctions = "results/{sample}/fetchdata/tumor_specific_sj.tsv",
#        tumor_associated_junctions = "results/{sample}/fetchdata/tumor_associated_sj.tsv",
#        expressed_in_healthy_junctions = "results/{sample}/fetchdata/expressed_in_healthy_sj.tsv"
#    params:
#        exe = workflow.source_path('../scripts/cancer_junction_selection.R'),
#        output_dir = lambda wildcards, output: os.path.dirname(output.tumor_specific_junctions),
#        healthy_tissue_sample_rate = 0.01
#    log:  "results/{sample}/log/cancer_junction_selection.log"
#    threads: 1
#    conda: '../envs/R.yaml'
#    shell:
#        'Rscript --vanilla {params.exe} '
#        '--sj {input.sj} '
#        '--filter_list {input.filter_list} '
#        '--filter_annotation {input.filter_list_annotation} '
#        '--healty_rate {params.healthy_tissue_sample_rate} '
#        '--output_dir {params.output_dir} 2>&1 | tee {log} '

rule find_cts_matching_wt:
    input:
        requantified_sj = "results/{sample}/fetchdata/requantified_sj.tsv",
        transcripts = os.path.join(config['index_dir'], 'ref_cdna.fa')
    output:
        blasted_junctions = "results/{sample}/fetchdata/requantified_blast_sj.tsv"
    params:
        tmp_dir = lambda wildcards, output: os.path.join(os.path.dirname(output[0]), 'blast_tmp'),
        exe = workflow.source_path('../scripts/identify_fp_requantification.py')
    threads: 8
    conda: '../envs/blast.yaml'
    log: "results/{sample}/log/blast.log"
    shell:
        'mkdir -p {params.tmp_dir}; '
        'python {params.exe} ' 
        '--input_splice {input.requantified_sj} '
        '--output {output} '
        '--transcriptome {input.transcripts} '
        '--temp_dir {params.tmp_dir} '
        '--threads {threads} 2>&1 | tee {log}'