# Input

TronMake RNA-splicing takes tumor RNA sequencing data as input. The input files are provided as tab separated value sheets (TSV) described in the following section. The pipeline supports multiple input modes / formats.

## Paired-end FASTQ

The table with paired end FASTQ files expects three tab-separated columns **without** a header.

| Sample name          | FASTQ 1                         | FASTQ 2                        |
|----------------------|---------------------------------|--------------------------------|
| sample_1             | /path/to/sample_1.1.fastq       |    /path/to/sample_1.2.fastq   |
| sample_2             | /path/to/sample_2.1.fastq       |    /path/to/sample_2.2.fastq   |

### Interleaved paired-end FASTQ

When using interleaved FASTQ input, the table expects two tab-separated columns **without** a header.
Make sure to set `interleaved_input: true` in the config file to trigger derinterleaving of the fastq file.

| sample   | interleaved fastq          |
|:--------:|:--------------------------:|
| sample_1 | /path/to/sample_1.fastq.gz |
| sample_2 | /path/to/sample_2.fastq.gz |


#### (u)BAM

When using (u)BAM files as input the table expects two tab-separated columns **without** a header.
Make sure to set `BAM_input: true` in the config file to trigger BAM to FASTQ conversion.

| Sample name   | BAM                                                   |
|:--------:|:-------------------------------------------------------:|
| sample_1 | /path/to/sample_1.bam                                   |
| sample_2 | /path/to/sample_2.bam,/path/to/sample_2_2.bam           |
| sample_3 | /path/to/sample_3.bam,/path/to/sample_3_2.bam           |

#### SRA 

When using the SRA mode, the table expects a single column **without** a header.
Make sure to set `sra_mode: true` in the config file to trigger download of FASTQ files of SRA Run accessions.

|accession |
|:--------:|
|SRR6298258| 
|SRR6298259|
|SRR6298260|
|SRR6298261|
|SRR6298262|

**Note: The pipeline currently support only SRA Run accessions starting with `SRR`, `DRR` or `ERR`. Other accessions such as study, expirment and/or group are currently not supported**
