# Output

The output of the pipeline consists of two main files.

* `results/{sample}/fetchdata/sj_final.tsv`
* `results/{sample}/fetchdata/sj_final_neofox_annotation.tsv`

Within the files, each line describes a candidate junction transcript. The file `sj_final.tsv` contains all candidate junctions with annotated features and targeted re-quantification results. The file `sj_final_neofox_annotation.tsv` contains only junction candidates with a mutated peptide sequence in format
suitable for analysis with NeoFox.

## Column-description sj_final.tsv

* `junc_id`: Predicted splice junction in standardised format `<chr>:<start>-<end>:<strand>`
* `uniquely_mapping_reads`: Uniquely mapped reads support by STAR.
* `multi_mapping_reads`: Multimapping reads support by STAR.
* `is_canonical`:  Logical indicating if predicted splice junction is in the set of reference junctions ("known junctions").
* `intron_jaccard`: Splice junction usage as defined by FRASER.
* `psi5`: PSI value of the donor splice site as defined by FRASER.
* `psi3`: PSI value of the acceptor splice site as defined by FRASER.
* `encode_blacklist_classification`: Junction is located ENCODE problematic region. Always NA in this output.
* `ucsc_blacklist_classification`: Junction is located UCSC problematic region. Always NA in this output.
* `tx_id`
* `tx_lst`
* `putative_event_type`
* `gene_id`
* `hgnc`
* `exclude_gene`
* `tx_mod_id`
* `junc_pos_tx`: The junction position in the modified transcript sequence.
* `cts_seq`: The context sequence. Thre transcript sequence in a defined window arround the junction position.
* `cts_junc_pos`:  The junction position in the context sequence.
* `cts_size`
* `cts_id`
* `transcript_expression_tpm`
* `gene_expression_tpm`
* `splice_site_motif`
* `jCPM_uniquely_mapped`
* `jCPM_multi_mapped`
* `jCPM_total_mapped`
* `junc_interval_start`
* `junc_interval_end`
* `span_interval_start`
* `span_interval_end`
* `within_interval`
* `within_interval_left`
* `within_interval_right`
* `coverage_perc`
* `coverage_perc_left`
* `coverage_perc_right`
* `coverage_mean`
* `coverage_mean_left`
* `coverage_mean_right`
* `coverage_median`
* `coverage_median_left`
* `coverage_median_right`
* `interval`
* `interval_left`
* `interval_right`
* `protein`
* `protein_wt`
* `frame_shift`
* `cds_mod_id`
* `junc_pos_cds`
* `protein_junc_pos`
* `is_first_reading_frame`
* `normalized_cds_junc_pos`
* `normalized_protein_junc_pos`
* `junc_in_orf`
* `truncated_cds`
* `cds_description`
* `peptide_context_seq_raw`
* `peptide_context_junc_pos`
* `peptide_context`
