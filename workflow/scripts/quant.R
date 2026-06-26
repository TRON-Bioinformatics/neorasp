suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(argparse)
  library(splice2neo)
})

df <- readr::read_tsv(snakemake@input[["sj"]], show_col_types = FALSE)

df <-
  splice2neo::map_requant(
    path_to_easyquant_folder = snakemake@params[["requant_dir"]],
    junc_tib = df
  )

df %>% readr::write_tsv(snakemake@output[["requantified_sj"]])
