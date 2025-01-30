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
* `tx_id`: The Ensembl identifier of the transcript.
* `tx_lst`: 
* `putative_event_type`: Putative event that the junction causes in the transcript.
* `gene_id`: The Ensembl identifier of the gene.
* `hgnc`: The HGNC identifier of the gene.
* `exclude_gene`: Logical if junction is excluded because of gene.
* `tx_mod_id`: The modified transcript id.
* `junc_pos_tx`: The junction position in the modified transcript sequence.
* `cts_seq`: The context sequence. Thre transcript sequence in a defined window arround the junction position.
* `cts_junc_pos`:  The junction position in the context sequence.
* `cts_size`: The size of the context sequence.
* `cts_id`: A unique id for the context sequence as hash value using the XXH128 hash algorithm.
* `transcript_expression_tpm`: Salmon expression estimate of gene in TPM.
* `gene_expression_tpm`: Salmon expression estimate of gene in TPM.
* `splice_site_motif`: Di-nucleotide motif of splice junction.
* `jCPM_uniquely_mapped`: CPM normalized uniquely mapped reads (STAR).
* `jCPM_multi_mapped`: CPM normalized multimapping reads (STAR).
* `jCPM_total_mapped`: CPM normalized total read support (unique + multimapping).
* `junc_interval_start`: Targeted re-quantification reads overlapping the junction of interest.
* `junc_interval_end`: Always NA.
* `span_interval_start`: Read pairs spanning the junction of interest.
* `span_interval_end`: Always NA.
* `within_interval`
* `within_interval_left`: Number of reads that map to interval left of splice junction of interest. 
* `within_interval_right`: Number of reads that map to interval right of splice junction of interest. 
* `coverage_perc`: Always NA.
* `coverage_perc_left`: Percentage of the interval left to splice junction that is covered by reads. 
* `coverage_perc_right`: Percentage of the interval right to splice junction that is covered by reads. 
* `coverage_mean`: Always NA.
* `coverage_mean_left`: Mean number of reads covering a position in the interval left to the splice junction of interest. 
* `coverage_mean_right`: Mean number of reads covering a position in the interval right to splice junction of interest. 
* `coverage_median`: Median number of reads covering a position in the interval left ti splice junction of interest. 
* `coverage_median_left`: Median number of reads covering a position in the interval left to splice junction of interest. 
* `coverage_median_right`: Median number of reads covering a position in the interval left to splice junction of interest.
* `interval`: Always NA
* `interval_left`: The interval of the context sequence left of the junction position the context sequence. 
* `interval_right`: The interval of the context sequence right of the junction position the context sequence. 
* `protein`:  The full protein sequence of the translated modified CDS.
* `protein_wt`:  The full protein sequence of the wild-type (reference) transcript.
* `frame_shift`
* `cds_mod_id`
* `junc_pos_cds`
* `protein_junc_pos`
* `is_first_reading_frame`
* `normalized_cds_junc_pos`
* `normalized_protein_junc_pos`
* `junc_in_orf`
* `truncated_cds`:  Indicator whether the mutated protein is a truncated from of the WT. If TRUE, `peptide_context` = NA.
* `cds_description`
* `peptide_context_seq_raw`
* `peptide_context_junc_pos`
* `peptide_context`
