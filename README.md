# NeoRasp

<!-- badges: start -->

[![Snakemake](https://img.shields.io/badge/snakemake-9.13.7-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
![Python](https://img.shields.io/badge/python-3670A0?style=flat-square&logo=python&logoColor=ffdd54)
[![R](https://img.shields.io/badge/R-276DC3?style=flat-square&logo=r&logoColor=white)](https://img.shields.io/badge/R-276DC3?style=flat-square&logo=r&logoColor=white)
[![CI](https://github.com/TRON-Private/tronmake-rna-splicing/actions/workflows/ci.yml/badge.svg)](https://github.com/TRON-Private/tronmake-rna-splicing/actions/workflows/ci.yml/badge.svg)
[![MIT](https://img.shields.io/badge/MIT-green?style=flat)](https://img.shields.io/badge/MIT-green?style=flat)

<!-- badges: end -->

**(Neo)antigens from (R)n(a)-(sp)licing**

______________________________________________________________________

<img align="right" width="150" height="150" src="https://github.com/user-attachments/assets/7f359faf-8a20-42c5-949d-4b2a70137d0f">  

Full documentation: [https://tron.pages.gitlab.rlp.net/tronmake-rna-splicing](https://probable-barnacle-2nrmw6z.pages.github.io/)

NeoRasp is an end-to-end workflow to identify non-canonical tumor-specific splice junction from RNA-seq.
The workflow implements a sensitive alignment based splice junction detection and targeted re-quantification
of candidate transcript variants. In our bioinformatics pipeline, [SnakeMake](https://snakemake.readthedocs.io/en/stable/) is employed as the primary workflow manager to orchestrate various steps.

______________________________________________________________________

## Workflow

- Input:

  - A table with paired-end FASTQ data for tumor samples.
  - A reference genome library. See [TronMake Genome Lib Builder](https://gitlab.rlp.net/tron/tronmake-genome-lib-builder)

- Process:

  1. Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp))

  1. Detection and metric calculation: ([`STAR`](https://github.com/alexdobin/STAR) -> [`fraser`](https://github.com/gagneurlab/FRASER))

  1. Expression quantification ([`Salmon`](https://salmon.readthedocs.io/en/latest/))

  1. Filtering:

     - Filtering based on junction expression.
     - Removing canonical junctions from GENCODE and healthy long read studies.
     - Removing junctions located in:
       - Problematic regions (low mappability).
       - IG-, TCR-, BCR- and HLA-regions.

  1. Targeted re-quantification of splice junction candidates ([`easyquant`](https://github.com/TRON-Bioinformatics/easyquant))

  1. Peptide annotation for [`NeoFox`](https://github.com/TRON-Bioinformatics/neofox)

## Dependencies

- python (>=3.10)
- snakemake (==9.13.7)
- conda (>=24.9)
- apptainer (>=1.3.4)

## Installation

### Download

```
git clone https://gitlab.rlp.net/tron/tronmake-rna-splicing

```

### Create conda environment

```
cd neorasp
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

- Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635.Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.
- Ines F. Scheller, Karoline Lutz, Christian Mertes, Vicente A. Yépez, Julien Gagneur medRxiv 2023.03.31.23287997; doi: https://doi.org/10.1101/2023.03.31.23287997
- Franziska Lang, Patrick Sorn, Martin Suchan, Alina Henrich, Christian Albrecht, Nina Köhl, Aline Beicht, Pablo Riesgo-Ferreiro, Christoph Holtsträter, Barbara Schrörs, David Weber, Martin Löwer, Ugur Sahin, Jonas Ibn-Salem, Prediction of tumor-specific splicing from somatic mutations as a source of neoantigen candidates, Bioinformatics Advances, Volume 4, Issue 1, 2024, vbae080, https://doi.org/10.1093/bioadv/vbae080
