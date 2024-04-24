library(tidyverse)
library(FRASER)
library(splice2neo)

saveFdsAsCountTable <- function(fds, min_read=5) {
    # junction info
    junction_dt <- as.data.table(granges(rowRanges(fds)))
    junction_dt$width <- NULL # remove uncessary columns
    colnames(junction_dt) <- c("Chromosome", "Start", "End", "Strand") # rename to fulfill the format https://gitlab.cmm.in.tum.de/gagneurlab/count_table
    junction_dt$Start <- as.integer(junction_dt$Start - 1) # change 1 based to 0 based. see https://www.biostars.org/p/84686/
    junction_dt$Strand <- gsub("[*]",".", junction_dt$Strand) # replace * with .
    # junction count
    count_dt <- tibble(raw_count = as.vector(counts(fds)))
    psi5 <- tibble(psi5 = as.vector(fds@assays@data$psi5))
    psi3 <- tibble(psi3 = as.vector(fds@assays@data$psi3))
    count_csv <- tibble(junction_dt) %>%
        bind_cols(count_dt) %>%
        bind_cols(psi5) %>%
        bind_cols(psi3) %>%
        filter(raw_count >= min_read)
    return(count_csv)
}

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

sample <- snakemake@wildcards[['sample']]
bam <- snakemake@input[['bam']]

sample_sheet <-
    tibble(sampleID = sample, bamFile = bam, gene = NA, pairedEnd = TRUE )
# Generate FRASER setting object
settings <- 
    FraserDataSet(
        colData = as.data.table(sample_sheet), 
         workingDir=snakemake@params[['working_dir']])

settings@bamParam@mapqFilter <- as.integer(snakemake@params[['mapq_filter']])
# Strand specific analysis 
strandSpecific(settings) <- as.integer(1)
# Count reads in BAM files with minimal read filtering
fds <- 
    countRNAData(settings,
        recount = TRUE,
        minExpressionInOneSample = snakemake@params[['min_read']],
        keepNonStandardChromosomes = FALSE)
# Calculate PSI values for splice junctions
fds <- calculatePSIValues(fds)
# Extract 5 and 3 PSI value for each junction
fds <- 
    saveFdsAsCountTable(fds, min_read = snakemake@params[['min_read']]) %>%
    dplyr::mutate(junc_id = splice2neo::generate_junction_id(Chromosome, Start, End, Strand))
# Write output
fds %>% write_tsv(snakemake@output[['psi_table']])

