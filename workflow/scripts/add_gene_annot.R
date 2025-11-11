suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(splice2neo)
  library(GenomicFeatures)
  library(furrr)
  library(purrr)
})



options(future.fork.enable = TRUE)
plan(multicore, workers = as.integer(snakemake@threads))

transcripts <- base::readRDS(snakemake@input[['transcripts']])
# Get transcripts with a defined CDS
cds <- base::readRDS(snakemake@input[['cds']])
cds_transcripts <- names(cds)

df <- readr::read_tsv(snakemake@input[['parsed_sj']], show_col_types = FALSE)

# Read transcript to gene mapping
tx2gene <- readr::read_tsv(snakemake@input[['tx2gene']], show_col_types = FALSE) %>%
  dplyr::rename(tx_id = TXNAME, gene_id = GENEID)

# Read gene to HGNC mapping
gene2hgnc <- readr::read_tsv(snakemake@input[['gene2hgnc']], show_col_types = FALSE) %>%
  dplyr::rename(hgnc = gene_symbol) %>%
  dplyr::distinct()

# Read RMSK annotation
rmsk <- base::readRDS(snakemake@input[['rmsk']])

# Annotate with possible transcripts
df <- df %>%
    splice2neo::add_tx(transcripts = transcripts) %>%
    splice2neo::annotate_potential_jet(., rmsk)

df_without_tx <- df %>%
    filter(is.na(tx_id)) %>%
    dplyr::mutate(hgnc="", gene_id="")
  
df_nc_transcripts <- df %>%
  dplyr::filter(!is.na(tx_id)) %>%
  dplyr::filter(!tx_id %in% cds_transcripts)

# Apply choose_tx in parallel to fix issue in splice2neo with large dataframes.
# 1. Remove junctions not overlapping any annotated transcript
# 2. Filter for junctions overlapping a transcript with defined CDS
df <- df %>%
  dplyr::filter(!is.na(tx_id)) %>%
  dplyr::filter(tx_id %in% cds_transcripts) %>%
  base::split(.$junc_id) %>% 
  furrr::future_map(~splice2neo::choose_tx(.x))
df <- dplyr::bind_rows(discard(df, ~nrow(.x) == 0))

# Select likely transcripts and annotate with ENSEMBL gene id and HGNC symbol
df <- df %>%
  dplyr::left_join(tx2gene) %>%
  dplyr::left_join(gene2hgnc)

df_without_tx %>%
  readr::write_tsv(snakemake@output[['annotated_sj_problematic']])

df_nc_transcripts %>%
  readr::write_tsv(snakemake@output[['annotated_sj_non_coding']])

df %>% readr::write_tsv(snakemake@output[['annotated_sj']])
