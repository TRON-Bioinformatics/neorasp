suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(magrittr)
  library(dplyr)
  library(argparse)
  library(splice2neo)
})

df <- readr::read_tsv(snakemake@input[["sj"]], show_col_types = FALSE)

# Filter only for consensus calls and remove intergenic calls
df %>%
  splice2neo::transform_for_requant() %>%
  readr::write_tsv(snakemake@output[["easyquant_table"]])

genes <- df %>%
  dplyr::pull(gene_id) %>%
  unique()

genes <- paste(genes, collapse = "|")

genes %>% readr::write_lines(snakemake@output[["genes_of_interest"]])
