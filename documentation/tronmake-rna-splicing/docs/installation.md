# Installation

This sections describes the process to install the workflow on a linux machine.
Please make sure the following dependencies are installed on your system. snakemake
and python are installed by the conda environment shipped with the repository.
Apptainer is only required when running the workflow with container support.

## Download pipeline repository

```
git clone https://gitlab.rlp.net/tron/tronmake-rna-splicing
```

## Create conda environment

```
cd tronmake-rna-splicing
conda env create -f environment.yaml --prefix conda_env/
conda activate conda_env
```

# System dependencies

 - python (>=3.10)
 - snakemake (==8.24.1)
 - conda (>=24.9)
 - apptainer (>=1.3.4)

# Pipeline dependencies

SnakeMake comes with integrated package management to retrieve and install all software
required to run the pipeline. The following table gives an overview which conda envrionments or
Docker containers are used by individual steps in the pipeline.

{{ read_table('../assets/software.tsv', sep = '\t') }}
