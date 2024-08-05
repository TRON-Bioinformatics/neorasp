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
parser$add_argument('--output', help= 'Output file')

xargs<- parser$parse_args()


df <- readr::read_tsv(xargs$sj)

tx <- readr::read_tsv(xargs$txp) %>%
  dplyr::rename(tx_id = Name, transcript_expression_tpm=TPM) %>%
  dplyr::select(tx_id, transcript_expression_tpm)

df <- df %>%
  dplyr::left_join(tx)

df %>% 
  readr::write_tsv(xargs$output)
