suppressMessages({
  options(stringsAsFactors = F)
  library(readr)
  library(dplyr)
  library(magrittr)
  library(argparse)
  library(splice2neo)
  library(rtracklayer)
  library(Biostrings)
  library(stringr)
  library(purrr)
})

# Read snakemake input/parameters
df <- readr::read_tsv(snakemake@input[['annotated_sj']], show_col_types = FALSE)
cds <- base::readRDS(snakemake@input[['cds']])
bsg <- rtracklayer::TwoBitFile(snakemake@input[['genome']])
peptide_flank_size <- as.integer(snakemake@params[["peptide_flank_size"]])
# Add splice junctions with peptides
peptide_annot <- df %>% 
  dplyr::select(junc_id, tx_id, cts_seq, cts_junc_pos, cts_size, cts_id) %>% 
  dplyr::distinct()

alt_peptides <- peptide_annot %>%
  splice2neo::add_peptide(cds = cds, flanking_size = peptide_flank_size, bsg = bsg)

# Calculate protein id for peptide FASTA
alt_peptides <- alt_peptides %>%
  dplyr::mutate(
    protein_id = 
      purrr::map2_chr(protein, protein_junc_pos,
      ~ {
        if (is.na(.x) || is.na(.y)) {
          NA_character_
        } else {
          rlang::hash(list(.x, .y))
        }
      })
  )

# Remove junctions without peptide altering effect
df <- df %>% 
  dplyr::left_join(alt_peptides) %>%
  dplyr::filter(cds_description == "mutated cds")

# Generate table for NeoFox annotation
dat_for_neofox <- df %>%
  dplyr::filter(!is.na(peptide_context) & nchar(peptide_context) > 7) %>%
  dplyr::filter(cds_description == "mutated cds") %>%
  dplyr::mutate(
    mutatedXmer = peptide_context,
    wildTypeXmer = NA,
    patientIdentifier = snakemake@wildcards[['sample']],
    rnaExpression = transcript_expression_tpm,
    rnaVariantAlleleFrequency = intron_jaccard,
    gene = hgnc,
  ) %>% dplyr::select(patientIdentifier, junc_id, tx_id, mutatedXmer, wildTypeXmer, rnaExpression, rnaVariantAlleleFrequency, gene) %>%
  dplyr::distinct()

# Write output files
df %>% readr::write_tsv(snakemake@output[['peptide_junc']])
dat_for_neofox %>% readr::write_tsv(snakemake@output[['neofox_annotation']])

# Assemble the FASTA header for Ligandomics analysis
df <- df %>%
  dplyr::mutate(
    fasta_header = paste0(
      "db_rna|",
      stringr::str_c("splice_", protein_id),
      "|",
      hgnc,
      " ",
      stringr::str_c("splice_", protein_id),
      " OS=Homo sapiens OX=9606 GN=",
      hgnc
    )
  )

# Remove junctions without peptide
df_fasta <- df %>%
  dplyr::select(fasta_header, protein, protein_junc_pos) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(protein) & !is.na(protein_junc_pos))

# Write FASTA output
peptides <- AAStringSet(df_fasta$protein)
names(peptides) <- df_fasta$fasta_header
writeXStringSet(peptides, snakemake@output[["peptide_fasta"]])
