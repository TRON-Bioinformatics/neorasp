# NeoRasp

<!-- badges: start -->

[![Release](https://gitlab.rlp.net/tron/tronmake-rna-splicing/-/badges/release.svg)](https://gitlab.rlp.net/tron/tronmake-rna-splicing/-/releases)
[![Snakemake](https://img.shields.io/badge/snakemake-9.1.3-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![pipeline status](https://gitlab.rlp.net/tron/tronmake-rna-splicing/badges/develop/pipeline.svg)](https://gitlab.rlp.net/tron/tronmake-rna-splicing/commits/master)

<!-- badges: end -->

NeoRasp in an end-to-end workflow to identify non-canonical tumor-specific splice junction from RNA-seq.
The workflow implements sensitive alignment-based detection of splice junctions and targeted re-quantification of candidate context sequences. 
In our bioinformatics pipeline, [SnakeMake](https://snakemake.readthedocs.io/en/stable/) is employed as the primary workflow manager to orchestrate various steps. 

The main steps include:  

  * Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp))  
  * Alignment of reads with STAR  ([`STAR`](https://github.com/alexdobin/STAR))  
  * Splice junction usage quantification ([`fraser`](https://github.com/gagneurlab/FRASER)) 
  * Gene and transcript expression quantification: ([`Salmon`](https://salmon.readthedocs.io/en/latest/))
  * Filtering of unwanted junctions.
  * Targeted re-quantification of splice junction candidates ([`easyquant`](https://github.com/TRON-Bioinformatics/easyquant))
  * Peptide annotation for neoantigen feature annotation ([`NeoFox`](https://github.com/TRON-Bioinformatics/neofox))


![NeoRasp DAG](assets/pipeline_rulegraph.svg){: style="height:700px;width:700px"}

		

