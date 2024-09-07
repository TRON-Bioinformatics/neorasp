# TronMake RNA-splice

<!-- badges: start -->

[![Release](https://gitlab.rlp.net/tron/tronmake-rna-splicing/-/badges/release.svg)](https://gitlab.rlp.net/tron/tronmake-rna-splicing/-/releases)
[![Snakemake](https://img.shields.io/badge/snakemake-8.16.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![pipeline status](https://gitlab.rlp.net/tron/tronmake-rna-splicing/badges/develop/pipeline.svg)](https://gitlab.rlp.net/tron/tronmake-rna-splicing/commits/master)

<!-- badges: end -->


The TronMake RNA-splice pipeline is a workflow to detect tumor-specific splice junctions in RNA-seq.

The workflow implements a sensitive alignment based splice junction detection and targeted re-quantification
of candidate transcript variants.

This pipeline aligns paired FASTQ files with STAR in 2-pass mode. It then calculates PSI and Intron jaccard values with FRASER.
Next, junctions are filtered for novelty. To exclude false-positve junction calls, we exclude alignments from problematic 
regions and genes. Remaining junction candidates are then re-quantified using easyquant.


## Components

1. Adapter and quality trimming ([`fastp`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))

2. Detection and PSI calculation: [`STAR`](https://github.com/alexdobin/STAR) -> [`fraser`](https://github.com/deweylab/RSEM) -> **SJ QUANTIFICATION**

3. Expression quantification: [`Salmon`](https://salmon.readthedocs.io/en/latest/) -> **GENE and TRANSCRIPT QUANTIFICATION**

4. Filtering: 
    * Removing canonical junctions from GENCODE, ENCODE, CHESS3 and GTEx long read GTFs
    * Removing junctions located in UCSC & ENCODE problematic regions.
    * Removing junctions from IG, TCR, BCR and HLA regions.

5. Targeted re-quantification of splice junctions ([`easyquant`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)).

6. Identification of false-positive re-quantification results.




## How to run it

Download the project and run as follows:
```
snakemake -p --use-conda --conda-prefix /path/to/shared_conda_envs \
    --executor local,slurm \
    --rerun-triggers mtime 
    --directory /path/to/output
    --configfile /path/to/config
```

### Input tables

The table with FASTQ files expects three tab-separated columns without a header

| Sample name          | FASTQ 1                      | FASTQ 2                  |
|----------------------|---------------------------------|------------------------------|
| sample_1             | /path/to/sample_1.1.fastq      |    /path/to/sample_1.2.fastq   |
| sample_2             | /path/to/sample_2.1.fastq      |    /path/to/sample_2.2.fastq   |

When using uBAM files as input the table expects two tab-separated columns without a header.

| Sample name   | BAM                                                   |
|:--------:|:-------------------------------------------------------:|
| sample_1 | /path/to/sample_1.bam                                   |
| sample_2 | /path/to/sample_2.bam,/path/to/sample_2_2.bam           |
| sample_3 | /path/to/sample_3.bam,/path/to/sample_3_2.bam           |


### Reference genome

The reference genome has to be provided as pre-built genome library. The TronMake RNA-expression
pipeline provides a script to create the required reference genome + tools indices.


## References


* Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635.Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.  
* Ines F. Scheller, Karoline Lutz, Christian Mertes, Vicente A. Yépez, Julien Gagneur medRxiv 2023.03.31.23287997; doi: https://doi.org/10.1101/2023.03.31.23287997   
* Franziska Lang, Patrick Sorn, Martin Suchan, Alina Henrich, Christian Albrecht, Nina Köhl, Aline Beicht, Pablo Riesgo-Ferreiro, Christoph Holtsträter, Barbara Schrörs, David Weber, Martin Löwer, Ugur Sahin, Jonas Ibn-Salem, Prediction of tumor-specific splicing from somatic mutations as a source of neoantigen candidates, Bioinformatics Advances, Volume 4, Issue 1, 2024, vbae080, https://doi.org/10.1093/bioadv/vbae080