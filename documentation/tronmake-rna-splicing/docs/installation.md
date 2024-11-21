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

