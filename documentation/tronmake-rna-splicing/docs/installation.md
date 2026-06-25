# Installation

This section describes how to install the workflow on Linux. The provided pixi environment installs Snakemake and Python. Apptainer is required for container support.

## Download pipeline repository

```
git clone https://gitlab.rlp.net/tron/tronmake-rna-splicing
```

## Install depedencies

```
pixi shell
```

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
pixi run test-local
```

## Test Structure

Tests are organized into different categories:

- **CI tests** (`--tag ci`): Run in GitHub Actions CI/CD pipeline
- **Local integration tests** (`--tag localintegrationtest`): Run on HPC systems with Apptainer
- **Issue-specific tests** (`tests/test_issue/`): Tests for specific bug fixes and features

See `tests/test_issue/README.md` for details on the issue-based testing convention.
