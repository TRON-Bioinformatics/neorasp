suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
})

df <- readr::read_tsv(snakemake@input[['annotated_sj']], show_col_types = FALSE)

tx <- readr::read_tsv(snakemake@input[['transcript_expression']], show_col_types = FALSE) %>%
  dplyr::rename(tx_id = Name, transcript_expression_tpm=TPM) %>%
  dplyr::select(tx_id, transcript_expression_tpm)

gene <- readr::read_tsv(snakemake@input[['gene_expression']], show_col_types = FALSE) %>%
  dplyr::rename(gene_id = Name, gene_expression_tpm=TPM) %>%
  dplyr::select(gene_id, gene_expression_tpm)

jx <- readr::read_tsv(snakemake@input[['junction_expression']], show_col_types = FALSE)

df <- df %>%
  dplyr::left_join(tx) %>%
  dplyr::left_join(gene) %>%
  dplyr::left_join(jx)

df %>% 
  readr::write_tsv(snakemake@output[['sj_expression']])
