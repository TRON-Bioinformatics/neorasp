suppressMessages({
  options(stringsAsFactors = F)
  library(FRASER)
  library(tidyverse)
  library(argparse)
  library(splice2neo)
})

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
    jaccard <- tibble(intron_jaccard = as.vector(fds@assays@data$jaccard))
    count_csv <- tibble(junction_dt) %>%
        bind_cols(count_dt) %>%
        bind_cols(psi5) %>%
        bind_cols(psi3) %>%
        bind_cols(jaccard) %>%
        dplyr::filter(raw_count >= min_read)
    return(count_csv)
}

parser <- ArgumentParser(description='Run FRASER PSI')

parser$add_argument('--bam', '-b', help='BAM file')
parser$add_argument('--output_table', '-o', help= 'Output table after conversion')
parser$add_argument('--threads', '-t', help= 'Number of threads')
parser$add_argument('--strandedness', '-s', help= '0(None), 1(FR), 2(True-Seq/RF)')
parser$add_argument('--min_expression', help= 'Minimum number of reads to call expression')
parser$add_argument('--mapq', help= 'Discard reads with smaller MAPQ for PSI calculation')



xargs<- parser$parse_args()


register(MulticoreParam(xargs$threads))

sample_sheet <-
    tibble(sampleID = "sample1", bamFile = xargs$bam, gene = NA, pairedEnd = TRUE )
# Generate FRASER setting object
settings <- 
    FraserDataSet(
        colData = as.data.table(sample_sheet), 
         workingDir=dirname(xargs$output_table))

settings@bamParam@mapqFilter <- as.integer(xargs$mapq)
# Strand specific analysis 
strandSpecific(settings) <- as.integer(xargs$strandedness)
# Count reads in BAM files with minimal read filtering
fds <- 
    countRNAData(
        settings,
        recount = TRUE,
        minExpressionInOneSample = as.integer(xargs$min_expression),
        keepNonStandardChromosomes = FALSE,
        NcpuPerSample = as.integer(xargs$threads)
    )
# Calculate PSI values for splice junctions
fds <- calculatePSIValues(fds)
# Extract 5 and 3 PSI value for each junction
fds <- 
    saveFdsAsCountTable(fds, min_read = as.integer(xargs$min_expression)) %>%
    dplyr::mutate(End = End + 1) %>%
    dplyr::mutate(junc_id = splice2neo::generate_junction_id(Chromosome, Start, End, Strand))
# Write output
fds %>% write_tsv(xargs$output_table)
