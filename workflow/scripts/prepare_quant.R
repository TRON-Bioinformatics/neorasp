suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(magrittr)
  library(dplyr)
  library(argparse)
  library(splice2neo)
  
})

parser <- ArgumentParser(description='Transform SJ sequences for easyquant')

parser$add_argument('--sj', help='SJ')
parser$add_argument('--output', help= 'EasyQuant input table')
parser$add_argument('--output_genes', help= 'Genes of interest for requant')

xargs<- parser$parse_args()


df <- readr::read_tsv(xargs$sj)

# Filter only for consensus calls and remove intergenic calls
df %>%
  splice2neo::transform_for_requant() %>% 
  readr::write_tsv(xargs$output)

genes <- df %>% dplyr::pull(gene) %>% unique()

genes <- paste(genes, collapse="|")

genes %>% readr::write_lines(xargs$output_genes)