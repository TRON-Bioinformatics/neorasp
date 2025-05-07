suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
})

df <- readr::read_tsv(snakemake@input[['annotated_sj']], show_col_types = FALSE)

df %>%
  dplyr::filter(intron_jaccard < snakemake@params[['min_junction_usage']]) %>%
  readr::write_tsv(snakemake@output[['sj_low_expression']])

df %>%
  dplyr::filter(intron_jaccard >= snakemake@params[['min_junction_usage']]) %>%
  readr::write_tsv(snakemake@output[['sj_expression']])
