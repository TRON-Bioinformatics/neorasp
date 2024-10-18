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

1. Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp))

2. Detection and PSI calculation: [`STAR`](https://github.com/alexdobin/STAR) -> [`fraser`](https://github.com/deweylab/RSEM) -> **SJ QUANTIFICATION**

3. Expression quantification: [`Salmon`](https://salmon.readthedocs.io/en/latest/) -> **GENE and TRANSCRIPT QUANTIFICATION**

4. Filtering: 
    * Filtering based on junction expression (either raw count or CPM based) 
    * Removing canonical junctions from GENCODE and healthy long read studies.
    * Removing junctions located in UCSC & ENCODE problematic regions.
    * Removing junctions from IG, TCR, BCR and HLA regions.

5. Targeted re-quantification of splice junction candidates ([`easyquant`](https://github.com/TRON-Bioinformatics/easyquant)).

6. Identification of false-positive re-quantification results.

7. Peptide annotation for NeoFox


## Installation

### Download

```
git clone https://gitlab.rlp.net/tron/tronmake-rna-splicing
```

### Create conda environment

```
cd tronmake-rna-splicing
conda env create -f environment.yaml --prefix conda_env/
conda activate conda_env
```

### Run the pipeline

```
snakemake \
    -p --use-conda --conda-prefix /path/to/shared_conda_envs \
    --executor local,slurm \
    --rerun-triggers mtime \
    --directory /path/to/output \
    --configfile /path/to/config \
    --config sample_sheet=/path/to/input.tsv \
    --software-deployment-method conda
```

* `--conda-prefix`: Where should the conda environments be stored
* `--config`: (Optional). Allows to overwrite settings from `--configfile`
* `--directory`: Specifies where the results are stored.
* `--configfile`: The path to the config file that contains e.g. the paths to the genome indices.
* `--software-deployment-method`: Conda and Apptainer are supported. For Apptainer support please make sure to bind input directories (Files + genome indices) into the container.


### Input tables

#### FASTQ

The table with normal paired end FASTQ files expects three tab-separated columns **without** a header

| Sample name          | FASTQ 1                      | FASTQ 2                  |
|----------------------|---------------------------------|------------------------------|
| sample_1             | /path/to/sample_1.1.fastq      |    /path/to/sample_1.2.fastq   |
| sample_2             | /path/to/sample_2.1.fastq      |    /path/to/sample_2.2.fastq   |


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

When using the SRA mode, the table expects two tab-separated columns **without** a header.
Make sure to set `sra_mode: true` in the config file to trigger download of FASTQ files of SRA Run accessions.

|accession |
|:--------:|
|SRR6298258| 
|SRR6298259|
|SRR6298260|
|SRR6298261|
|SRR6298262|

**Note: The pipeline currently support only SRA Run accessions starting with `SRR`. Other accessions such as study, expirment and/or group are not supported**



### Reference genome

The reference genome has to be provided as pre-built genome library. The (TronMake Genome Lib Builder)[https://gitlab.rlp.net/tron/tronmake-genome-lib-builder]
pipeline provide the required genome  and tool indices.


## Authors & Acknowledgements 

The TronMake RNA splicing pipeline was originally developed by Johannes Hausmann at [TRON - Translational Oncology at the Medical Center of the Johannes Gutenberg University Mainz gGmbH (non-profit)](https://tron-mainz.de/).

Maintenance is now lead by Johannes Hausmann. 

Main developers: 

- [Johannes Hausmann](mailto:johannes.hausmann@tron-mainz.de)   

Contributers:

- None / NA

## References

* Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635.Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.  
* Ines F. Scheller, Karoline Lutz, Christian Mertes, Vicente A. Yépez, Julien Gagneur medRxiv 2023.03.31.23287997; doi: https://doi.org/10.1101/2023.03.31.23287997   
* Franziska Lang, Patrick Sorn, Martin Suchan, Alina Henrich, Christian Albrecht, Nina Köhl, Aline Beicht, Pablo Riesgo-Ferreiro, Christoph Holtsträter, Barbara Schrörs, David Weber, Martin Löwer, Ugur Sahin, Jonas Ibn-Salem, Prediction of tumor-specific splicing from somatic mutations as a source of neoantigen candidates, Bioinformatics Advances, Volume 4, Issue 1, 2024, vbae080, https://doi.org/10.1093/bioadv/vbae080