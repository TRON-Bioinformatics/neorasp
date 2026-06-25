## Output Files

The pipeline generates the following main output files for each sample:

- `results/{sample}/star/{sample}_Aligned.sortedByCoord.out.cram` - Aligned reads in CRAM format
- `results/{sample}/fetchdata/sj_final.tsv` - Final junction candidates with annotations
- `results/{sample}/fetchdata/sj_final_neofox_annotation.tsv` - NeoFox-compatible annotation
- `results/{sample}/fetchdata/sj_final_peptides.fasta` - Altered peptide sequences
- `results/report/multiqc.html` - Quality control report

## Main Result Files

The output of the pipeline consists of three main files.

- `results/{sample}/fetchdata/sj_final.tsv`
- `results/{sample}/fetchdata/sj_final_neofox_annotation.tsv`
- `results/{sample}/fetchdata/sj_final_peptides.fasta`

Within the files, each line describes a candidate junction transcript. The file `sj_final.tsv` contains all candidate junctions with annotated features and targeted re-quantification results. The file `sj_final_neofox_annotation.tsv` contains only junction candidates with a mutated peptide sequence in format
suitable for analysis with NeoFox. `sj_final_peptides.fasta` containts the altered full-length protein sequence for MassSpec analysis. 

## Column-descriptions

### sj_final.tsv

- `junc_id`: Predicted splice junction in standardised format `<chr>:<start>-<end>:<strand>`
- `uniquely_mapping_reads`: Uniquely mapped reads support by STAR.
- `multi_mapping_reads`: Multimapping reads support by STAR.
- `is_canonical`: Splice junction is in the set of reference junctions from the genome library (boolean).
- `intron_jaccard`: Splice junction usage as defined by FRASER.
- `psi5`: PSI value of the donor splice site as defined by FRASER.
- `psi3`: PSI value of the acceptor splice site as defined by FRASER.
- `encode_blacklist_classification`: Junction is located ENCODE problematic region. Always NA in this output.
- `ucsc_blacklist_classification`: Junction is located UCSC problematic region. Always NA in this output.
- `tx_id`: The Ensembl identifier of the transcript.
- `putative_event_type`: Putative event that the junction causes in the transcript.
- `gene_id`: The Ensembl identifier of the gene.
- `hgnc`: The HGNC identifier of the gene.
- `junc_pos_tx`: The junction position in the modified transcript sequence.
- `cts_seq`: The context sequence. The transcript sequence in a defined window around the junction position.
- `cts_junc_pos`: The junction position in the context sequence.
- `cts_size`: The size of the context sequence.
- `cts_id`: A unique id for the context sequence as hash value using the XXH128 hash algorithm.
- `transcript_expression_tpm`: Salmon expression estimate of transcript in TPM.
- `gene_expression_tpm`: Salmon expression estimate of gene in TPM.
- `splice_site_motif`: Di-nucleotide motif of splice junction.
- `jCPM_uniquely_mapped`: CPM normalized uniquely mapped reads (STAR).
- `jCPM_multi_mapped`: CPM normalized multimapping reads (STAR).
- `jCPM_total_mapped`: CPM normalized total read support (unique + multimapping).
- `junction_reads`: Targeted re-quantification reads overlapping the junction of interest.
- `spanning_reads`: Read pairs spanning the junction of interest.
- `within_interval_left`: Number of reads that map to interval left of splice junction of interest.
- `within_interval_right`: Number of reads that map to interval right of splice junction of interest.
- `coverage_perc_left`: Percentage of the interval left to splice junction that is covered by reads.
- `coverage_perc_right`: Percentage of the interval right to splice junction that is covered by reads.
- `coverage_mean_left`: Mean number of reads covering a position in the interval left to the splice junction.
- `coverage_mean_right`: Mean number of reads covering a position in the interval right to splice junction.
- `coverage_median`: Median number of reads covering a position in the interval left ti splice junction.
- `coverage_median_left`: Median number of reads covering a position in the interval left to splice junction.
- `coverage_median_right`: Median number of reads covering a position in the interval left to splice junction.
- `interval_left`: The interval of the context sequence left of the junction position the context sequence.
- `interval_right`: The interval of the context sequence right of the junction position the context sequence.
- `protein`: The full protein sequence of the translated modified CDS.
- `protein_wt`: The full protein sequence of the wild-type (reference) transcript.
- `frame_shift`: Splice junction leads to frameshift in transcript (boolean).
- `junc_pos_cds`: The junction position in the modified CDS sequence.
- `protein_junc_pos`: The position of the junction in the protein sequence.
- `is_first_reading_frame`: Modified CDS sequence is translated into protein sequence using the 1st reading frame (boolean).
- `normalized_cds_junc_pos`: The normalized position of the junction in the modified CDS sequence to the left junction side.
- `normalized_protein_junc_pos`: The normalized position of the junction in the protein sequence to the left junction side.
- `junc_in_orf`: Splice junction is located in an open reading frame (boolean).
- `truncated_cds`: Mutated protein is a truncated from of the WT (boolean). If TRUE, `peptide_context` = NA.
- `cds_description`: Descriptor of of the mutated gene product. Can be one of c("mutated cds", "truncated cds", "no mutated gene product", "no wt cds", "not in ORF")
- `peptide_context_seq_raw`: The peptide sequence around the junction including stop codons.
- `peptide_context_junc_pos`: The junction position relative to the `peptide_context` sequence.
- `peptide_context`: The peptide sequence around the junction truncated after stop codons.

### sj_final_neofox_annotation.tsv

- `patientIdentifier`: Patient identifier for NeoFox input. Must match to `identifier` in NeoFox patient data sheet.
- `mutatedXmer`: The peptide sequence around the junction.
- `wildTypeXmer`: Peptide sequence of wild-type. Currently always NA.
- `rnaExpression`: Transcript expression in TPM.
- `rnaVariantAlleleFrequency`: Matches to `intron_jaccard`.
- `gene`: The HGNC identifier of the gene.

### sj_final_peptides.fasta

This fasta file containts the altered full-length peptide sequence for MassSpec analysis. The header of each entry has the following format.

- `>db_rna|splice_<protein_id>|<hgnc> splice_<protein_id> OS=Homo sapiens OX=9606 GN=<hgnc>`

- `protein_id`: Unique id for the protein sequence as hash value using the XXH128 hash algorithm.

- `hgnc`: The HGNC identifier of the gene.

An example heaader from the test data:

- `>db_rna|splice_8760f9325ca45b33c224399f28e374ea|DNAI7 splice_8760f9325ca45b33c224399f28e374ea OS=Homo sapiens OX=9606 GN=DNAI7`
