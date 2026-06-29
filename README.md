<p align="center">
    <img src="https://github.com/user-attachments/assets/7f359faf-8a20-42c5-949d-4b2a70137d0f" alt="logo" width="150" height="150">
</p>



# NeoRasp

<!-- badges: start -->

[![Snakemake](https://img.shields.io/badge/snakemake-9.13.7-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
![Python](https://img.shields.io/badge/python-3670A0?style=flat-square&logo=python&logoColor=ffdd54)
[![R](https://img.shields.io/badge/R-276DC3?style=flat-square&logo=r&logoColor=white)](https://img.shields.io/badge/R-276DC3?style=flat-square&logo=r&logoColor=white)
[![CI](https://github.com/TRON-Private/tronmake-rna-splicing/actions/workflows/ci.yml/badge.svg)](https://github.com/TRON-Private/tronmake-rna-splicing/actions/workflows/ci.yml/badge.svg)
[![MIT](https://img.shields.io/badge/MIT-green?style=flat)](https://img.shields.io/badge/MIT-green?style=flat)

<!-- badges: end -->

**(Neo)antigens from (R)n(a)-(sp)licing**

Full documentation: [https://tron-bioinformatics.github.io/neorasp/](https://tron-bioinformatics.github.io/neorasp/)
______________________________________________________________________

**NeoRasp** is an end-to-end workflow for identifying non-canonical, tumor-specific intra-gene splice junctions from RNA-seq data. It employs sensitive alignment-based detection of splice junctions followed by targeted re-quantification of candidate context sequences. The pipeline is orchestrated using [Snakemake](https://snakemake.readthedocs.io/en/stable/).


## Workflow

- Process:

  1. Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp))

  2. Detection and metric calculation: ([`STAR`](https://github.com/alexdobin/STAR) -> [`fraser`](https://github.com/gagneurlab/FRASER))

  3. Expression quantification ([`Salmon`](https://salmon.readthedocs.io/en/latest/))

  4. Filtering:

     - Filtering based on junction expression.
     - Removing canonical junctions from GENCODE and healthy long read studies.
     - Removing junctions located in:
       - Problematic regions (low mappability).
       - IG-, TCR-, BCR- and HLA-regions.

  5. Targeted re-quantification of splice junction candidates ([`easyquant`](https://github.com/TRON-Bioinformatics/easyquant))

  6. Peptide annotation for [`NeoFox`](https://github.com/TRON-Bioinformatics/neofox)

## Dependencies for installation

- pixi
- apptainer (>=1.3.4)

## Installation

### Download

```
git clone https://github.com/TRON-Bioinformatics/neorasp
```

### Install dependencies with pixi

```
pixi shell
```

## Usage

To run NeoRasp, adapt and execute the following command:

```
snakemake -s workflow/Snakefile \
	--directory </path/to/output/directory> \
    --latency-wait 60 \
    --software-deployment-method apptainer \
	[--configfile <path/to/config/file>] \
	[--profile </path/to/cluster/profile/>]
```

## Input

NeoRasp requires user-provided input:

  1. A table with paired-end FASTQ data for tumor samples. See [*Input section*](https://tron-bioinformatics.github.io/neorasp/input) of documentation.
  2. A reference genome library. See [OBLX genome library](https://github.com/TRON-Bioinformatics/oblx)

## Output

The output of the pipeline is written to the directory specified with
`--directory`. Descriptions of the output files are documented in
[*Output section*](https://tron-bioinformatics.github.io/neorasp/output) of the documentation.

## Contribution

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) and
[*developer_guide section*](https://tron-bioinformatics.github.io/neorasp/developer_guide/) of the documentation.


## Authors & Acknowledgements

The NeoRasp pipeline was originally developed in the Computational Genomics group at [TRON - Translational Oncology at the Medical Center of the Johannes Gutenberg University Mainz gGmbH (non-profit)](https://tron-mainz.de/).

Maintenance is currently led by Johannes Hausmann.

🛠️ Main developers:

- [Johannes Hausmann, TRON gGmbH](https://github.com/johausmann)  

✨ Contributors and 🐞 bug hunter:

- [Luis Kress, TRON gGmbH](https://github.com/LKress)
- [Franziska Lang, TRON gGmbH](https://github.com/franla23)
- [Jonas Ibn-Salem, TRON gGmbH](https://github.com/ibn-salem)

## References

- Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635.Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.
- Ines F. Scheller, Karoline Lutz, Christian Mertes, Vicente A. Yépez, Julien Gagneur medRxiv 2023.03.31.23287997; doi: https://doi.org/10.1101/2023.03.31.23287997
- Franziska Lang, Patrick Sorn, Martin Suchan, Alina Henrich, Christian Albrecht, Nina Köhl, Aline Beicht, Pablo Riesgo-Ferreiro, Christoph Holtsträter, Barbara Schrörs, David Weber, Martin Löwer, Ugur Sahin, Jonas Ibn-Salem, Prediction of tumor-specific splicing from somatic mutations as a source of neoantigen candidates, Bioinformatics Advances, Volume 4, Issue 1, 2024, vbae080, https://doi.org/10.1093/bioadv/vbae080
- Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Van Dyken, P. C., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Vieira, F. G., Meesters, C., Lee, S., Twardziok, S. O., Kanitz, A., VanCampen, J., Malladi, V., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2025). Sustainable data analysis with Snakemake. F1000Research, 10, 33. https://doi.org/10.12688/f1000research.29032.3