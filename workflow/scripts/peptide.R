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



df <- readr::read_tsv(snakemake@input[['sj']], show_col_types = FALSE)
cds <- base::readRDS(snakemake@input[['cds']])
bsg <- rtracklayer::TwoBitFile(snakemake@input[['genome']])

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

df %>% readr::write_tsv(snakemake@output[['junctions']])
dat_for_neofox %>% readr::write_tsv(snakemake@output[['neofox_annotation']])