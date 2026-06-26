# Cluster execution

## Custom Profile

For execution on clusters managed by workload managers like SLURM, a profile file tailored to the cluster's configuration is required. A generic example for SLURM is provided below:

```
executor: cluster-generic
cluster-generic-submit-cmd: sbatch -n {threads} --mem={resources.mem_mb} --job-name=smk-{rule} --parsable
cluster-generic-cancel-cmd: scancel
jobs: 100
```

The file has to be located in a directory called e.g. slurm and the file has to be named `config.v8+.yaml`.
To use the profile, just call the pipeline with the parameter `--profile <path/to/slurm_folder>`.
The jobs parameter in the profile file specifies the maximum number of jobs that Snakemake can submit in parallel.
When the pipeline is started within a SLURM job, no more than 1 core and 1 GB of memory is required for this initial
job as all tasks are submitted as individual jobs. However, we recommend to start the initial SnakeMake process
from an interactive shell session.

## SnakeMake executor plugin

SnakeMake provides multiple executor plugins to run the workflow on compute clusters. So far, we have only tested the slurm-executor-plugin.
