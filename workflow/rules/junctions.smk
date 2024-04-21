rule fraser:
    input:
        bam = get_aligner_sam,
        gtf = get_genome_annotation
    output:
        psi_table = ""
    script:
        '../scripts/fraser.R'

rule superintronic:
    input:
        bam = ,
        gtf = ,
    output:

rule splicemap:
    pass

rule regtools:
    pass

rule splice2neo:
    pass

rule requant:
    pass

rule annotate_paralogus:
    pass

saveFdsAsCountTable <- function(fds, min_read) {
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
