# Input

NeoRasp takes tumor RNA sequencing data as input. The input files are provided as tab separated value sheets (TSV) described in the following section. The pipeline supports multiple input modes / formats.

## Sequencing data

### Paired-end FASTQ

The table with paired end FASTQ files expects three tab-separated columns **without** a header. This
is the reommended input data type.

| Sample name | FASTQ 1                   | FASTQ 2                   |
| ----------- | ------------------------- | ------------------------- |
| sample_1    | /path/to/sample_1.1.fastq | /path/to/sample_1.2.fastq |
| sample_2    | /path/to/sample_2.1.fastq | /path/to/sample_2.2.fastq |

## Genome library

The workflow requires the genome annotation as well as indices for STAR and Salmon as genome library which are specified in the config file. Please use TronMake genome library for analysis.
