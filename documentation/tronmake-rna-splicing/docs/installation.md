# Installation

This sections describes the process to install the workflow on a linux machine.
Snakemake and python are installed by the conda environment shipped with the repository.
Apptainer is required to run the workflow with container support.

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

**Please make sure that the conda chanel priority is not set to strict.**
Channel priority is configured and set in `$HOME/.condarc`.

# System dependencies

Please make sure the following dependencies are installed on your system.

 - conda (>=24.9)
 - apptainer (>=1.3.4)

