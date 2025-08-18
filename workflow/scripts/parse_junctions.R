suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(argparse)
  library(splice2neo)
  library(tidyr)
  library(fs)
})

parse_star_junctions <- function(star_sj, fraser_sj, canonical_junctions, cpm_sj) {
  df_fraser <-
    readr::read_tsv(fraser_sj, 
                     col_types = cols(Start = col_integer(), End = col_integer(), Strand=col_character()),
                     show_col_types = FALSE) %>%
    dplyr::rename(junction_start = Start, 
                  junction_end = End,
                  strand = Strand,
                  chromosome = Chromosome,
                  number_supporting_reads = raw_count)

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

  canonical_juncs <- readr::read_tsv(canonical_junctions, show_col_types = FALSE)
  
  cpm_sj <- readr::read_tsv(cpm_sj, show_col_types = FALSE)

  df_star <- splice2neo::parse_star_sj(path = star_sj) %>%
    dplyr::mutate(is_canonical = junc_id %in% canonical_juncs$junc_id) %>%
    dplyr::left_join(cpm_sj) 

  filtered_junctions <- df_star %>%
    dplyr::filter(is_canonical) %>%
    dplyr::select(-c(Gene, class))

  df_star <- df_star %>%
    dplyr::filter(!is_canonical) %>% dplyr::select(-c(Gene, class))

  # Filter out junctions with a supporting read count less than the provided cutoff
  df_star <- df_star %>%
    dplyr::left_join(df_fraser %>% dplyr::select(junc_id, intron_jaccard, psi5, psi3) %>% distinct())

  return(list(novel_junctions = df_star, canonical_junctions = filtered_junctions)) 

}

junctions <- parse_star_junctions(
    star_sj = snakemake@input[['star_sj']],
    fraser_sj = snakemake@input[['fraser_psi']],
    canonical_junctions = snakemake@input[['canonical_junctions']],
    cpm_sj = snakemake@input[['star_cpm']]
)

junctions$novel_junctions %>%
  readr::write_tsv(snakemake@output[['parsed_sj']])
junctions$canonical_junctions %>%
  readr::write_tsv(snakemake@output[['removed_junction']])
