library(fraser)
library(tidyverse)

samples <- snakemake@wilcards[[1]]
bam <- snakemake@input[['bam']]

