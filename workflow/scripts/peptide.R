suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(argparse)
  library(splice2neo)
  library(rtracklayer)
})

parser <- ArgumentParser(description='Annotate junctions with peptide sequence')

parser$add_argument('--sj', help='Splice junctions')
parser$add_argument('--sample', help= 'Sample name')
parser$add_argument('--cds', help= 'Sample name')
parser$add_argument('--genome', help= 'Sample name')
parser$add_argument('--output', help= 'Output')
parser$add_argument('--output_neofox', help= 'Output for neofox')

xargs<- parser$parse_args()

df <- readr::read_tsv(xargs$sj, show_col_types = FALSE)
cds <- base::readRDS(xargs$cds)
bsg <- rtracklayer::TwoBitFile(xargs$genome)

peptide_annot <- df %>% 
  dplyr::select(junc_id, tx_id, cts_seq, cts_junc_pos, cts_size, cts_id) %>% 
  dplyr::distinct()

alt_peptides <- peptide_annot %>%
  splice2neo::add_peptide(cds = cds, flanking_size = 13, bsg = bsg)

df <- df %>% 
  dplyr::left_join(alt_peptides)

dat_for_neofox <- df %>%
  dplyr::filter(!is.na(peptide_context) & nchar(peptide_context) > 7) %>%
  dplyr::mutate(
    mutatedXmer = peptide_context,
    wildTypeXmer = NA,
    patientIdentifier = xargs$sample
  )

df %>% readr::write_tsv(xargs$output)
dat_for_neofox %>% readr::write_tsv(xargs$output_neofox)