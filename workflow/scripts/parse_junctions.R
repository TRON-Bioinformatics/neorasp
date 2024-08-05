suppressMessages({
  options(stringsAsFactors = F)
  library(tidyverse)
  library(argparse)
  library(splice2neo)
  library(GenomicFeatures)
})

parser <- ArgumentParser(description='Parse junctions and annotate with splice2neo')

parser$add_argument('--sj', help='STAR predicted splice junctions')
parser$add_argument('--fraser', help='FRASER PSI and jaccard values')
parser$add_argument('--read_support', help='Number of spliced alignments to keep a junction.')
parser$add_argument('--output', help= 'Output directory')

xargs<- parser$parse_args()

canonical_juncs <- readr::read_tsv(xargs$canonical_juncs)


df_fraser <-
    read_tsv(xargs$fraser,  col_types = cols(Start = col_character(), End = col_character())) %>%
    dplyr::rename(junction_start = Start, 
                  junction_end = End,
                  strand = Strand,
                  chromosome = Chromosome,
                  number_supporting_reads = raw_count) %>%
    dplyr::filter(chromosome %in% paste0('chr', c(1:22, 'X', 'Y')))
  
# Annotate junctins with plus and minus strand when unstranded data and let STAR decide the strand of the junctions based on XS tag *
df_fraser <- df_fraser %>%
  dplyr::select(junction_start, junction_end, strand, chromosome, intron_jaccard) %>%
  dplyr::mutate(
    junc_id = dplyr::case_when(
      strand == '*' ~ stringr::str_glue("{splice2neo::generate_junction_id(chromosome, start, end, '+')};{splice2neo::generate_junction_id(chromosome, start, end, '-')}"),
      TRUE ~ splice2neo::generate_junction_id(chromosome, start, end, strand)
    )
  ) %>% tidyr::separate_rows(junc_id, sep=";")

df_star <- splice2neo::parse_star_sj(path = xargs$canonical_juncs)

# Filter out junctions with a supporting read count less than the provided cutoff
star_sj <- star_sj %>%
  dplyr::mutate(number_supporting_reads = as.numeric(uniquely_mapping_reads) + as.numeric(multi_mapping_reads)) %>%
  dplyr::filter(number_supporting_reads >= as.numeric(xargs$read_support)) %>%
  dplyr::filter(chromosome %in% paste0('chr', c(1:22, 'X', 'Y'))) %>% 
  dplyr::left_join(df_fraser %>% dplyr::select(junc_id, intron_jaccard, psi5, psi3))




  
  