# TronMake Cancer RNA-splicing

<!-- badges: start -->

[![Release](https://gitlab.rlp.net/tron/tronmake-rna-splicing/-/badges/release.svg)](https://gitlab.rlp.net/tron/tronmake-rna-splicing/-/releases)
[![Snakemake](https://img.shields.io/badge/snakemake-8.24.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![pipeline status](https://gitlab.rlp.net/tron/tronmake-rna-splicing/badges/develop/pipeline.svg)](https://gitlab.rlp.net/tron/tronmake-rna-splicing/commits/master)

<!-- badges: end -->

TronMake Cancer RNA-splicing is a workflow to identify non-canonical tumor-specific splice junction from RNA-seq.
The workflow implements a sensitive alignment based splice junction detection and targeted re-quantification
of candidate transcript variants. In our bioinformatics pipeline, [SnakeMake](https://snakemake.readthedocs.io/en/stable/) is employed as the primary workflow manager to orchestrate various steps. 

The main steps comprise:  

  * Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp))  
  * Alignment of reads with STAR  ([`STAR`](https://github.com/alexdobin/STAR))  
  * Splice junction usage quantification ([`fraser`](https://github.com/gagneurlab/FRASER)) 
  * Gene and transcript expression quantification: ([`Salmon`](https://salmon.readthedocs.io/en/latest/))
  * Filtering of spurious junctions.
  * Targeted re-quantification of splice junction candidates ([`easyquant`](https://github.com/TRON-Bioinformatics/easyquant))
  * Peptide annotation for neoantigen feature annotation ([`NeoFox`](https://github.com/TRON-Bioinformatics/neofox))


## Components
|Tool | Link | License|
|:----:|:----:|:------:|
|fastp| https://github.com/OpenGene/fastp | MIT|
|STAR | https://github.com/alexdobin/STAR | MIT|
|samtools |https://github.com/samtools/samtools | MIT/Expat|
|fraser	| https://github.com/gagneurlab/FRASER | MIT|
|salmon	| https://salmon.readthedocs.io/en/latest/ | GPL-3.0|
|easyquant |https://github.com/TRON-Bioinformatics/easyquant | MIT|
|splice2neo | https://github.com/TRON-Bioinformatics/splice2neo | MIT|
		

