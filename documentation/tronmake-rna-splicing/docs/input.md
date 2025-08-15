# Input

NeoRasp takes tumor RNA sequencing data as input. The input files are provided as tab separated value sheets (TSV) described in the following section. The pipeline supports multiple input modes / formats.

## Sequencing data

### Paired-end FASTQ

The table with paired end FASTQ files expects three tab-separated columns **without** a header. This
is the reommended input data type.

| Sample name          | FASTQ 1                         | FASTQ 2                        |
|----------------------|---------------------------------|--------------------------------|
| sample_1             | /path/to/sample_1.1.fastq       |    /path/to/sample_1.2.fastq   |
| sample_2             | /path/to/sample_2.1.fastq       |    /path/to/sample_2.2.fastq   |

#### Interleaved paired-end FASTQ

When using interleaved FASTQ input, the table expects two tab-separated columns **without** a header.
Make sure to set `interleaved_input: true` in the config file to trigger derinterleaving of the fastq file.

| sample   | interleaved fastq          |
|:--------:|:--------------------------:|
| sample_1 | /path/to/sample_1.fastq.gz |
| sample_2 | /path/to/sample_2.fastq.gz |



## Genome library

The workflow requires the genome annotation as well as indices for STAR and Salmon as genome library which are specified in the config file. In future, the libraries of the TronMake Genome Lib Builder will be used. 
