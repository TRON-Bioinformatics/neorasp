suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(argparse)
  library(splice2neo)
  
})

parser <- ArgumentParser(description='Annotate junctions with easyquant requant results')

parser$add_argument('--sj', help='spladder')
parser$add_argument('--requant', help='requant')
parser$add_argument('--output', help= 'Output directory')

xargs<- parser$parse_args()


df <- readr::read_tsv(xargs$sj, show_col_types = FALSE)

df <-
    splice2neo::map_requant(
        path_to_easyquant_folder = xargs$requant,
        junc_tib = df
    )

df <- df %>%
  dplyr::mutate(junction_psi = intron_jaccard * gene_expression_tpm)

df %>% readr::write_tsv(xargs$output)
