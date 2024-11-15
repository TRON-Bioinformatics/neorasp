suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(argparse)
  library(splice2neo)
  library(GenomicFeatures)
  requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
})

gene_exclusion <- function(df, exclusion_pattern) {
  #stopifnot("hgnc" %in% colnames(df), "Column HGNC is missing")
  exclusion_pattern <- readr::read_tsv(exclusion_pattern, show_col_types = FALSE)
  
  df <- df %>%
    fuzzyjoin::regex_left_join(exclusion_pattern,
      by = c("hgnc" = "exclude_gene_pattern")
    ) %>%
    dplyr::mutate(exclude_gene =!is.na(exclude_gene_pattern))
  return(df)

}

transcripts <- base::readRDS(snakemake@input[['transcripts']])
bsg <- rtracklayer::TwoBitFile(snakemake@input[['genome']])

df <- readr::read_tsv(snakemake@input[['star_sj']], show_col_types = FALSE)

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
    dplyr::mutate(hgnc="", gene_id="", exclude_gene=FALSE)

df <- df %>%
  splice2neo::choose_tx()

# Select likely transcripts and annotate with ENSEMBL gene id and HGNC symbol
df <- df %>%
  dplyr::left_join(tx2gene) %>%
  dplyr::left_join(gene2hgnc)

df <- gene_exclusion(df, snakemake@input[['gene_exclusion']])

df %>% 
  dplyr::filter(exclude_gene | is.na(gene_id)) %>%
  dplyr::bind_rows(df_without_tx) %>%
  readr::write_tsv(nakemake@output[['annotated_sj_problematic']])

# Merge with gene and remove genes from exclusion pattern

df <- df %>%
  dplyr::filter(!exclude_gene) %>%
  dplyr::filter(!is.na(gene_id)) %>%
  splice2neo::add_context_seq(transcripts = transcripts, size = snakemake@params[['cts_size']], bsg = bsg)

df %>% readr::write_tsv(snakemake@output[['annotated_sj']])
