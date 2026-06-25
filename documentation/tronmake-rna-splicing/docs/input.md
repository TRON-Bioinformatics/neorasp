# Input

NeoRasp takes tumor RNA-seq data as input. The input files are provided as tab separated value sheets (TSV) described in the following section.

## Sequencing data

The table with paired end FASTQ files expects three tab-separated columns **without** a header. 

| Sample name          | FASTQ 1                         | FASTQ 2                        |
|----------------------|---------------------------------|--------------------------------|
| sample_1             | /path/to/sample_1.1.fastq       |    /path/to/sample_1.2.fastq   |
| sample_2             | /path/to/sample_2.1.fastq       |    /path/to/sample_2.2.fastq   |

### Technical replicates

NeoRasp supports the processing of technical replicates by providing them in the input table together. The table containing the replicate FASTQ files requires three tab-separated columns **without** a header. NeoRasp aligns the replicates together in one file and sets the SM tag in the BAM file to differentiate reads from the input files.

| Sample name          | FASTQ 1                         | FASTQ 2                        |
|----------------------|---------------------------------|--------------------------------|
| sample_1             | /path/to/sample_11.1.fastq,/path/to/sample_12.1.fastq       |    /path/to/sample_11.2.fastq,/path/to/sample_12.2.fastq   |
| sample_2             | /path/to/sample_21.1.fastq,/path/to/sample_22.1.fastq       |    /path/to/sample_21.2.fastq,/path/to/sample_22.2.fastq   |

## OBLX genome library

The workflow requires the genome annotation as well as indices for STAR and Salmon  which are specified in the config file. The [OBLX genome library](https://github.com/TRON-Bioinformatics/oblx) provides all files required to run the pipeline. If you want to use your own genome annotation adapt the [config](usage.md#config) file to your needs.
