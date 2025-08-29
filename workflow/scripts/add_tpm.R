suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(splice2neo)
})

df <- readr::read_tsv(snakemake@input[['annotated_sj']], show_col_types = FALSE)

tx <- readr::read_tsv(snakemake@input[['transcript_expression']], show_col_types = FALSE) %>%
  dplyr::rename(tx_id = Name, transcript_expression_tpm=TPM) %>%
  dplyr::select(tx_id, transcript_expression_tpm)

gene <- readr::read_tsv(snakemake@input[['gene_expression']], show_col_types = FALSE) %>%
  dplyr::rename(gene_id = Name, gene_expression_tpm=TPM) %>%
  dplyr::select(gene_id, gene_expression_tpm)

transfrags <- readr::read_tsv(snakemake@input[['transfrags_expression']], show_col_types = FALSE) %>%
  dplyr::mutate(junc_id = splice2neo::generate_junction_id(chrom, start - 1, end + 1, strand)) %>%
  dplyr::select(junc_id, tx_ids, stringtie_TPM) %>%
  dplyr::rename(stringtie_transfrags_ids = tx_ids, stringtie_transfrags_tpm = stringtie_TPM)

df <- df %>%
  dplyr::left_join(tx) %>%
  dplyr::left_join(gene) %>%
  dplyr::left_join(transfrags)

df %>% 
  readr::write_tsv(snakemake@output[['sj_expression']])
