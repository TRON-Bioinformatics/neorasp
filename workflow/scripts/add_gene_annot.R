suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(splice2neo)
  library(GenomicFeatures)
})

transcripts <- base::readRDS(snakemake@input[['transcripts']])

df <- readr::read_tsv(snakemake@input[['parsed_sj']], show_col_types = FALSE)

# Read transcript to gene mapping
tx2gene <- readr::read_tsv(snakemake@input[['tx2gene']], show_col_types = FALSE) %>%
  dplyr::rename(tx_id = TXNAME, gene_id = GENEID)

# Read gene to HGNC mapping
gene2hgnc <- readr::read_tsv(snakemake@input[['gene2hgnc']], show_col_types = FALSE) %>%
  dplyr::select(`Gene stable ID version`, `Gene name`) %>%
  dplyr::rename(gene_id = `Gene stable ID version`, hgnc = `Gene name`) %>%
  dplyr::distinct()

# Annotate with possible transcripts
df <- df %>%
    splice2neo::add_tx(transcripts = transcripts)

df_without_tx <- df %>%
    filter(is.na(tx_id)) %>%
    dplyr::mutate(hgnc="", gene_id="")

df <- df %>%
  splice2neo::choose_tx()

# Select likely transcripts and annotate with ENSEMBL gene id and HGNC symbol
df <- df %>%
  dplyr::left_join(tx2gene) %>%
  dplyr::left_join(gene2hgnc)

df_without_tx %>%
  readr::write_tsv(snakemake@output[['annotated_sj_problematic']])

df %>% readr::write_tsv(snakemake@output[['annotated_sj']])
