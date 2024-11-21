# TronMake cancer RNA-splicing

<!-- badges: start -->

[![Release](https://gitlab.rlp.net/tron/tronmake-rna-splicing/-/badges/release.svg)](https://gitlab.rlp.net/tron/tronmake-rna-splicing/-/releases)
[![Snakemake](https://img.shields.io/badge/snakemake-8.24.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![pipeline status](https://gitlab.rlp.net/tron/tronmake-rna-splicing/badges/develop/pipeline.svg)](https://gitlab.rlp.net/tron/tronmake-rna-splicing/commits/master)

<!-- badges: end -->

Full documentation: https://tron.pages.gitlab.rlp.net/tronmake-rna-splicing


The TronMake Cancer RNA-splicing pipeline is a workflow to detect tumor-specific splice junctions in tumor RNA.
The workflow implements a sensitive alignment based splice junction detection and targeted re-quantification
of candidate transcript variants. This pipeline aligns paired-end FASTQ files with STAR in 2-pass mode. It then calculates 
percent splice in (PSI) and Intron jaccard values with FRASER. Next, junctions are filtered for novelty. To exclude 
false-positve junction calls, we exclude alignments from problematic regions and genes. Lastly, remaining junction candidates will be re-quantified with easyquant.


## Workflow

- Input:
  - A table with paired-end FASTQ data for tumor samples
  - A reference genome lib

- Process:
    1. Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp))

    2. Detection and metric calculation: [`STAR`](https://github.com/alexdobin/STAR) -> [`fraser`](https://github.com/deweylab/RSEM) -> **SJ QUANTIFICATION**

    3. Expression quantification: [`Salmon`](https://salmon.readthedocs.io/en/latest/) -> **GENE and TRANSCRIPT QUANTIFICATION**

    4. Filtering: 
        * Filtering based on junction expression.
        * Removing canonical junctions from GENCODE and healthy long read studies.
        * Removing junctions located in:
            * Problematic regions (low mappability).
            * IG-, TCR-, BCR- and HLA-regions.
    
    5. Targeted re-quantification of splice junction candidates ([`easyquant`](https://github.com/TRON-Bioinformatics/easyquant)).

    6. Identification of false-positive re-quantification results.

    7. Peptide annotation for NeoFox

## Dependencies

 - python (>=3.10)
 - snakemake (==8.24.1)
 - conda (>=24.9)
 - apptainer (>=1.3.4)

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


## Authors & Acknowledgements 

The TronMake Cancer RNA splicing pipeline was originally developed by Johannes Hausmann at [TRON - Translational Oncology at the Medical Center of the Johannes Gutenberg University Mainz gGmbH (non-profit)](https://tron-mainz.de/).

Maintenance is now lead by Johannes Hausmann. 

Main developers: 

- [Johannes Hausmann](mailto:johannes.hausmann@tron-mainz.de)   

Contributers:

- None / NA / NULL

## References

* Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635.Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.  
* Ines F. Scheller, Karoline Lutz, Christian Mertes, Vicente A. Yépez, Julien Gagneur medRxiv 2023.03.31.23287997; doi: https://doi.org/10.1101/2023.03.31.23287997   
* Franziska Lang, Patrick Sorn, Martin Suchan, Alina Henrich, Christian Albrecht, Nina Köhl, Aline Beicht, Pablo Riesgo-Ferreiro, Christoph Holtsträter, Barbara Schrörs, David Weber, Martin Löwer, Ugur Sahin, Jonas Ibn-Salem, Prediction of tumor-specific splicing from somatic mutations as a source of neoantigen candidates, Bioinformatics Advances, Volume 4, Issue 1, 2024, vbae080, https://doi.org/10.1093/bioadv/vbae080