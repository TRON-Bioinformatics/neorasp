suppressMessages({
  options(stringsAsFactors = F)
  library(FRASER)
  library(magrittr)
  library(dplyr)
  library(tibble)
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

runFraserAnalysis <- function(sample_sheet, strandedness = 0, min_expression=5, mapq_filter=255, workingDir, threads) {

    settings <- 
        FraserDataSet(
            colData = as.data.table(sample_sheet), 
            workingDir=workingDir)

    settings@bamParam@mapqFilter <- mapq_filter

    strandSpecific(settings) <- as.integer(strandedness)

    fds <- countRNAData(
                settings,
                recount = TRUE,
                minExpressionInOneSample = as.integer(min_expression),
                keepNonStandardChromosomes = FALSE,
                NcpuPerSample = as.integer(threads)
            )
    # Calculate PSI values for splice junctions
    fds <- calculatePSIValues(fds)
    return(fds)
}




register(MulticoreParam(snakemake@threads))

sample_sheet <-
    tibble(sampleID = "sample1", bamFile = snakemake@input[['bam']], gene = NA, pairedEnd = TRUE)

fds <- runFraserAnalysis(sample_sheet, 
                         strandedness = 0,
                         min_expression = snakemake@params[['min_read']],
                         mapq_filter = snakemake@params[['mapq_filter']],
                         workingDir = dirname(snakemake@output[['psi_table']]),
                         threads = snakemake@threads)

fds <- saveFdsAsCountTable(fds, min_read = snakemake@params[['min_read']]) %>%
    dplyr::mutate(End = End + 1) %>%
    dplyr::mutate(junc_id = stringr::str_c(Chromosome, ":", Start,"-", End,":", Strand))
# Write output
fds %>% write.table(., file = snakemake@output[['psi_table']], sep="\t")
