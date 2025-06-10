suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(splice2neo)
  library(GenomicFeatures)
  library(purrr)
})

transcripts <- base::readRDS(snakemake@input[['transcripts']])
bsg <- rtracklayer::TwoBitFile(snakemake@input[['genome']])

df <- readr::read_tsv(snakemake@input[['parsed_sj']], show_col_types = FALSE)

df <- df %>%
  dplyr::filter(!exclude_gene) %>%
  dplyr::filter(!is.na(gene_id)) %>%
  splice2neo::add_context_seq(.,
      transcripts = transcripts,
      size = snakemake@params[['cts_size']],
      bsg = bsg)
  )

df %>% readr::write_tsv(snakemake@output[['annotated_sj']])
