suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(argparse)
  library(splice2neo)
  library(GenomicFeatures)
})

parser <- ArgumentParser(description='Parse junctions and annotate with splice2neo')

parser$add_argument('--sj', help='STAR predicted splice junctions')
parser$add_argument('--canonical_junctions', help='Canonical splice junctions to exclude from futher analysis')
parser$add_argument('--fraser', help='FRASER PSI and jaccard values')
parser$add_argument('--read_support', help='Number of spliced alignments to keep a junction.')
parser$add_argument('--output', help= 'Output directory')

xargs<- parser$parse_args()

df_fraser <-
    read_tsv(xargs$fraser,  col_types = cols(Start = col_integer(), End = col_integer(), Strand=col_character()), show_col_types = FALSE) %>%
    dplyr::rename(junction_start = Start, 
                  junction_end = End,
                  strand = Strand,
                  chromosome = Chromosome,
                  number_supporting_reads = raw_count) %>%
    dplyr::filter(chromosome %in% paste0('chr', c(1:22, 'X', 'Y')))

# Annotate junctins with plus and minus strand when unstranded data and let STAR decide the strand of the junctions based on XS tag *
df_fraser <- df_fraser %>%
  dplyr::select(junction_start, junction_end, strand, chromosome, intron_jaccard, psi5, psi3) %>%
  rowwise() %>%
  dplyr::mutate(
    junc_id = dplyr::case_when(
      strand == '.' ~ stringr::str_glue("{splice2neo::generate_junction_id(chromosome, junction_start, junction_end, '+')};{splice2neo::generate_junction_id(chromosome, junction_start, junction_end, '-')}"),
      strand == '+' ~ splice2neo::generate_junction_id(chromosome, junction_start, junction_end, '+'),
      strand == '-' ~ splice2neo::generate_junction_id(chromosome, junction_start, junction_end, '-')
    )
  ) %>% tidyr::separate_rows(junc_id, sep=";")

canonical_junctions <- readr::read_tsv(xargs$canonical_junctions, show_col_types = FALSE)

df_star <- splice2neo::parse_star_sj(path = xargs$sj) %>%
  dplyr::mutate(is_canonical = junc_id %in% canonical_junctions$junc_id) 

filtered_junctions <- df_star %>%
  dplyr::filter(is_canonical) %>%
  dplyr::select(-c(Gene, class)) %>% dplyr::left_join(canonical_junctions)

df_star <- df_star %>%
  dplyr::filter(!is_canonical) %>% dplyr::select(-c(Gene, class))

# Filter out junctions with a supporting read count less than the provided cutoff
df_star <- df_star %>%
  dplyr::filter(as.numeric(uniquely_mapping_reads) >= as.numeric(xargs$read_support)) %>%
  dplyr::filter(chromosome %in% paste0('chr', c(1:22, 'X', 'Y'))) %>% 
  dplyr::left_join(df_fraser %>% dplyr::select(junc_id, intron_jaccard, psi5, psi3) %>% distinct())

df_star %>% readr::write_tsv(fs::path(xargs$output, "parsed_sj.tsv.tmp"))
filtered_junctions %>% readr::write_tsv(fs::path(xargs$output, "sj_canonical.tsv"))



  
  