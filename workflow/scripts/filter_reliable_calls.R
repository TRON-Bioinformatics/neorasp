suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
})

df <- readr::read_tsv(snakemake@input[["annotated_sj"]], show_col_types = FALSE)

df <- df %>%
  dplyr::mutate(
    reliable_call = intron_jaccard >= snakemake@params[["min_junction_usage"]] &
      jCPM_uniquely_mapped >= snakemake@params[["min_junction_cpm"]]
  )

df %>%
  dplyr::filter(!reliable_call) %>%
  readr::write_tsv(snakemake@output[["sj_low_expression"]])

df %>%
  dplyr::filter(reliable_call) %>%
  readr::write_tsv(snakemake@output[["sj_expression"]])
