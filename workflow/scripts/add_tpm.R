suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(argparse)
  
})

parser <- ArgumentParser(description='Annotate junctions with transcript expression')

parser$add_argument('--sj', help='junctions')
parser$add_argument('--txp', help='transcript expression from salmon')
parser$add_argument('--gxp', help='gene expression from salmon')
parser$add_argument('--output', help= 'Output file')

xargs<- parser$parse_args()


df <- readr::read_tsv(xargs$sj, show_col_types = FALSE)

tx <- readr::read_tsv(xargs$txp, show_col_types = FALSE) %>%
  dplyr::rename(tx_id = Name, transcript_expression_tpm=TPM) %>%
  dplyr::select(tx_id, transcript_expression_tpm)

gene <- read_tsv(xargs$gxp, show_col_types = FALSE) %>%
  dplyr::rename(gene_id = Name, gene_expression_tpm=TPM) %>%
  dplyr::select(gene_id, gene_expression_tpm)


df <- df %>%
  dplyr::left_join(tx) %>%
  dplyr::left_join(gene)

df %>% 
  readr::write_tsv(xargs$output)
