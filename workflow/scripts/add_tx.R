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
#gene_exclusion <- function(df, exclusion_pattern) {
#  #stopifnot("hgnc" %in% colnames(df), "Column HGNC is missing")
#  exclusion_pattern <- readr::read_tsv(exclusion_pattern, show_col_types = FALSE)
#  
#  df <- df %>%
#    fuzzyjoin::regex_left_join(exclusion_pattern,
#      by = c("hgnc" = "exclude_gene_pattern")
#    ) %>%
#    dplyr::mutate(exclude_gene =!is.na(exclude_gene_pattern))
#  return(df)
#
#}

transcripts <- base::readRDS(snakemake@input[['transcripts']])
bsg <- rtracklayer::TwoBitFile(snakemake@input[['genome']])

df <- readr::read_tsv(snakemake@input[['parsed_sj']], show_col_types = FALSE)

df <- df %>%
  dplyr::filter(!exclude_gene) %>%
  dplyr::filter(!is.na(gene_id)) %>% 
  base::split(.$junc_id) %>% 
  furrr::future_map(
    ~splice2neo::add_context_seq(.x, 
      transcripts = transcripts,
      size = snakemake@params[['cts_size']],
      bsg = bsg)
  )
df <- dplyr::bind_rows(discard(df, ~nrow(.x) == 0))

df %>% readr::write_tsv(snakemake@output[['annotated_sj']])
