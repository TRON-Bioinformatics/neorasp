# Input

NeoRasp takes tumor RNA-seq data as input. The input files are provided as tab separated value sheets (TSV) described in the following section..

## Sequencing data

### Paired-end FASTQ

The table with paired end FASTQ files expects three tab-separated columns **without** a header. 

| Sample name          | FASTQ 1                         | FASTQ 2                        |
|----------------------|---------------------------------|--------------------------------|
| sample_1             | /path/to/sample_1.1.fastq       |    /path/to/sample_1.2.fastq   |
| sample_2             | /path/to/sample_2.1.fastq       |    /path/to/sample_2.2.fastq   |

### Paired-end FASTQ (technical replicates)

NeoRasp supports the processing of technical replicates by providing them in the input table together. The table containing the replicate FASTQ files requires three tab-separated columns **without** a header. NeoRasp aligns the replicates together in one file and sets the SM tag in the BAM file to differentiate reads from the input files.

| Sample name          | FASTQ 1                         | FASTQ 2                        |
|----------------------|---------------------------------|--------------------------------|
| sample_1             | /path/to/sample_11.1.fastq,/path/to/sample_12.1.fastq       |    /path/to/sample_11.2.fastq,/path/to/sample_12.2.fastq   |
| sample_2             | /path/to/sample_21.1.fastq,/path/to/sample_22.1.fastq       |    /path/to/sample_21.2.fastq,/path/to/sample_22.2.fastq   |

## Genome library

The workflow requires the genome annotation as well as indices for STAR and Salmon as genome library which are specified in the config file. Please use TronMake genome library for analysis.
