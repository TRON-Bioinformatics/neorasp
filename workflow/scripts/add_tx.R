suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(splice2neo)
  library(GenomicFeatures)
  library(purrr)
  library(fs)
  library(withr)
  library(stringr)
})

output_dir <- dirname(snakemake@output[["annotated_sj"]])

tmp_genome <- tempfile(
  pattern = stringr::str_glue("genome_copy_{snakemake@wildcards[['chunkID']]}_"),
  tmpdir = output_dir,
  fileext = ".2bit"
)

tmp_genome_final <- tempfile(
  pattern = stringr::str_glue("genome_copy_{snakemake@wildcards[['chunkID']]}_"),
  tmpdir = output_dir,
  fileext = ".2bit"
)

tmp_transcripts <- tempfile(
  pattern = stringr::str_glue("transcript_copy_{snakemake@wildcards[['chunkID']]}_"),
  tmpdir = output_dir,
  fileext = ".RDS"
)

# Copy to tmp
# Perform atomar NFS operation to ensure objects dont use same cache
fs::file_copy(snakemake@input[["genome"]], tmp_genome, overwrite = TRUE)
fs::file_move(tmp_genome, tmp_genome_final)

fs::file_copy(snakemake@input[["transcripts"]], tmp_transcripts, overwrite = TRUE)

defer({
  if (file_exists(tmp_genome_final)) file_delete(tmp_genome_final)
  if (file_exists(tmp_transcripts)) file_delete(tmp_transcripts)
})

bsg <- rtracklayer::TwoBitFile(tmp_genome_final)
transcripts <- base::readRDS(tmp_transcripts)

df <- readr::read_tsv(snakemake@input[["parsed_sj"]], show_col_types = FALSE)

df <- df %>%
  dplyr::filter(!exclude_gene) %>%
  dplyr::filter(!is.na(gene_id)) %>%
  splice2neo::add_context_seq(.,
    transcripts = transcripts,
    size = snakemake@params[["cts_size"]],
    bsg = bsg
  )

df %>% readr::write_tsv(snakemake@output[["annotated_sj"]])
