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
cd NeoRasp
conda env create -f environment.yaml --prefix conda_env/
conda activate conda_env
```

**Please make sure that the conda chanel priority is not set to strict. Channel priority is configured and set in `$HOME/.condarc`.**

# System dependencies

Please make sure the following dependencies are installed on your system.

 - conda (>=24.9)
 - apptainer (>=1.3.4)

# Testing the Pipeline

NeoRasp includes comprehensive integration tests to verify the pipeline functionality.

## Running Local Integration Tests

For local testing on HPC systems with Apptainer support, use the provided Makefile:

```bash
# Set the path to your Apptainer library
export APPTAINER_HPC=/path/to/apptainer/library

# Run all local integration tests
make localintegrationtest
```

The Makefile runs pytest with the `localintegrationtest` tag, which executes tests designed for local environments with proper Apptainer configurations.

## Test Structure

Tests are organized into different categories:

* **CI tests** (`--tag ci`): Run in GitHub Actions CI/CD pipeline
* **Local integration tests** (`--tag localintegrationtest`): Run on HPC systems with Apptainer
* **Issue-specific tests** (`tests/test_issue/`): Tests for specific bug fixes and features

See `tests/test_issue/README.md` for details on the issue-based testing convention.

