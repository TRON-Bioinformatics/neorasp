suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(magrittr)
  library(dplyr)
  library(argparse)
  library(rtracklayer)
  
})

filter_genomic_regions <- function(df, encode_blacklist, ucsc_unusal_regions){
  
  # build GRanges of junction
  #jx <- df %>% dplyr::select(junc_id) %>% distinct() %>% pull()
  jx <- splice2neo::junc_to_gr(df$junc_id)
  
  encode_blacklist <- rtracklayer::import.bed(encode_blacklist)
  # Find encode ranges that contain the junction range
  hits <- GenomicRanges::findOverlaps(jx, encode_blacklist, type = "within")
  jx_idx <- hits %>% S4Vectors::queryHits()
  encode_blacklist_idx <- hits %>% S4Vectors::subjectHits()
  # build data.frame to associate junc_id to transcripts
  junc_to_encode <- tibble::tibble(
    junc_id = df$junc_id[jx_idx],
    encode_blacklist_classification = encode_blacklist[encode_blacklist_idx]$name
  )

  ucsc_blacklist <- rtracklayer::import.bed(ucsc_unusal_regions)
  seqlevels(ucsc_blacklist, pruning.mode="coarse") <- seqlevelsInUse(jx)
  hits <- GenomicRanges::findOverlaps(jx, ucsc_blacklist, type = "within", ignore.strand=TRUE)
  jx_idx <- hits %>% S4Vectors::queryHits()
  ucsc_blacklist_idx <- hits %>% S4Vectors::subjectHits()
  
  # build data.frame to associate junc_id to transcripts
  junc_to_ucsc <- tibble::tibble(
    junc_id = df$junc_id[jx_idx],
    ucsc_blacklist_classification = ucsc_blacklist[ucsc_blacklist_idx]$name
  )
  
  
  # join the original input data.frame with the association of junction to transcripts
  out_df <- df %>%
    dplyr::left_join(junc_to_encode, by = "junc_id", relationship = "many-to-many") %>%
    dplyr::left_join(junc_to_ucsc, by = "junc_id", relationship = "many-to-many")
  
  return(out_df)
}

df <- readr::read_tsv(snakemake@input[['parsed_sj']], show_col_types = FALSE)

df <- df %>% 
  filter_genomic_regions(., 
    snakemake@input[['encode_regions']], snakemake@input[['ucsc_regions']])

df_failed <- df %>%
  dplyr::filter(!is.na(encode_blacklist_classification) | !is.na(ucsc_blacklist_classification))

df %>%
  dplyr::filter(is.na(encode_blacklist_classification) & is.na(ucsc_blacklist_classification)) %>%
  readr::write_tsv(snakemake@output[['parsed_sj']])

df_failed %>%
  readr::write_tsv(snakemake@output[['failed_sj']])
