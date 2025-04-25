# NeoRasp

<!-- badges: start -->

[![Release](https://gitlab.rlp.net/tron/tronmake-rna-splicing/-/badges/release.svg)](https://gitlab.rlp.net/tron/tronmake-rna-splicing/-/releases)
[![Snakemake](https://img.shields.io/badge/snakemake-9.1.3-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![pipeline status](https://gitlab.rlp.net/tron/tronmake-rna-splicing/badges/dev/pipeline.svg)](https://gitlab.rlp.net/tron/tronmake-rna-splicing/commits/main)

<!-- badges: end -->

**(Neo)antigens from (R)n(a)-(sp)licing**


Full documentation: https://tron.pages.gitlab.rlp.net/tronmake-rna-splicing


NeoRasp is an end-to-end workflow to identify non-canonical tumor-specific splice junction from RNA-seq.
The workflow implements a sensitive alignment based splice junction detection and targeted re-quantification
of candidate transcript variants. In our bioinformatics pipeline, [SnakeMake](https://snakemake.readthedocs.io/en/stable/) is employed as the primary workflow manager to orchestrate various steps. 

## Workflow

- Input:
  - A table with paired-end FASTQ data for tumor samples.
  - A reference genome library. See [TronMake Genome Lib Builder](https://gitlab.rlp.net/tron/tronmake-genome-lib-builder)

- Process:
    1. Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp))

    2. Detection and metric calculation: ([`STAR`](https://github.com/alexdobin/STAR) -> [`fraser`](https://github.com/deweylab/RSEM))

    3. Expression quantification ([`Salmon`](https://salmon.readthedocs.io/en/latest/))

    4. Filtering: 
        * Filtering based on junction expression.
        * Removing canonical junctions from GENCODE and healthy long read studies.
        * Removing junctions located in:
            * Problematic regions (low mappability).
            * IG-, TCR-, BCR- and HLA-regions.
    
    5. Targeted re-quantification of splice junction candidates ([`easyquant`](https://github.com/TRON-Bioinformatics/easyquant))

    6. Peptide annotation for [`NeoFox`](https://github.com/TRON-Bioinformatics/neofox)

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
cd NeoRasp
conda env create -f environment.yaml --prefix conda_env/
conda activate conda_env
```


## Authors & Acknowledgements 

The NeoRasp pipeline was originally developed by Johannes Hausmann at [TRON - Translational Oncology at the Medical Center of the Johannes Gutenberg University Mainz gGmbH (non-profit)](https://tron-mainz.de/).

Maintenance is now lead by Johannes Hausmann. 

Main developers: 

- [Johannes Hausmann](mailto:johannes.hausmann@tron-mainz.de)   

Contributers:

- Luis Kress, TRON gGmbH
- Franziska Lang, TRON gGmbH
- Jonas Ibn-Salem, TRON gGmbH

## References

* Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635.Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.  
* Ines F. Scheller, Karoline Lutz, Christian Mertes, Vicente A. Yépez, Julien Gagneur medRxiv 2023.03.31.23287997; doi: https://doi.org/10.1101/2023.03.31.23287997   
* Franziska Lang, Patrick Sorn, Martin Suchan, Alina Henrich, Christian Albrecht, Nina Köhl, Aline Beicht, Pablo Riesgo-Ferreiro, Christoph Holtsträter, Barbara Schrörs, David Weber, Martin Löwer, Ugur Sahin, Jonas Ibn-Salem, Prediction of tumor-specific splicing from somatic mutations as a source of neoantigen candidates, Bioinformatics Advances, Volume 4, Issue 1, 2024, vbae080, https://doi.org/10.1093/bioadv/vbae080