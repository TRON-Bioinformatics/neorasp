suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(argparse)
  library(splice2neo)
  library(rtracklayer)
  library(Biostrings)
  library(stringr)
  library(purrr)
})

df <- readr::read_tsv(snakemake@input[['sj']], show_col_types = FALSE)
cds <- base::readRDS(snakemake@input[['cds']])
bsg <- rtracklayer::TwoBitFile(snakemake@input[['genome']])

peptide_annot <- df %>% 
  dplyr::select(junc_id, tx_id, cts_seq, cts_junc_pos, cts_size, cts_id) %>% 
  dplyr::distinct()

alt_peptides <- peptide_annot %>%
  splice2neo::add_peptide(cds = cds, flanking_size = 13, bsg = bsg)

alt_peptides <- alt_peptides %>%
  dplyr::mutate(
    protein_id = 
      purrr::map2_chr(protein, protein_junc_pos,
      ~ {
        if (is.na(.x) || is.na(.y)) {
          NA_character_
        } else {
          rlang::hash(list(.x, .y))
        }
      })
  )

df <- df %>% 
  dplyr::left_join(alt_peptides) %>%
  dplyr::select(
    -tx_lst,
    -exclude_gene,
    -tx_mod_id,
    -junc_interval_end,
    -span_interval_end,
    -within_interval,
    -coverage_perc,
    -coverage_mean,
    -coverage_median,
    -interval,
    -cds_mod_id, 
  ) %>% dplyr::rename(junction_reads = junc_interval_start,
                      spanning_reads = span_interval_start)

dat_for_neofox <- df %>%
  dplyr::filter(!is.na(peptide_context) & nchar(peptide_context) > 7) %>%
  dplyr::mutate(
    mutatedXmer = peptide_context,
    wildTypeXmer = NA,
    patientIdentifier = snakemake@wildcards[['sample']],
    rnaExpression = transcript_expression_tpm,
    rnaVariantAlleleFrequency = intron_jaccard,
    gene = hgnc,
  ) %>% dplyr::select(patientIdentifier, junc_id, tx_id, mutatedXmer, wildTypeXmer, rnaExpression, rnaVariantAlleleFrequency, gene) %>%
  dplyr::distinct()

df %>% readr::write_tsv(snakemake@output[['junctions']])
dat_for_neofox %>% readr::write_tsv(snakemake@output[['neofox_annotation']])

df <- df %>%
  dplyr::mutate(
    fasta_header = paste0(
      "db_rna|",
      stringr::str_c("splice_", protein_id),
      "|",
      hgnc,
      " ",
      stringr::str_c("splice_", protein_id),
      " OS=Homo sapiens OX=9606 GN=",
      hgnc
    )
  )

df_fasta <- df %>%
  dplyr::select(fasta_header, protein) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(protein))

peptides <- AAStringSet(df_fasta$protein)
names(peptides) <- df_fasta$fasta_header
writeXStringSet(peptides, snakemake@output[["peptide_fasta"]])
