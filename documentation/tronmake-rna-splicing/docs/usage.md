# Usage

> **Information: If possible, use absolute paths when specifying directories or files.**  
> **Recommendation: If possible, configure the workflow with a yaml configfile and not on the command line.**

```
snakemake \
    --config sample_sheet=</path/to/samples.tsv> index_dir=</path/to/genome_lib> \
    --directory <output_directory> \
    --software-deployment-method conda \
    --configfile <path to my config> \
    --workflow-profile workflow/profiles/default \
    --default-resources "tmpdir='<output_directory>'" \
    [--conda-prefix </path/to/conda/env/location>] \
    [--profile </path/to/profile/directory>]
    [--retries <n>]
```

* `config`: After `--config` different parameters can be specified (**NOTE**: the parameters after config are separated by space and do not start with `--`. **Moreover they can be omitted if a configfile is used**):
    * `sample_sheet=`: Define the path to the sample sheet.
    * `index_dir=`: Path to the genome library and genome indices.
* `directory`: Specifies where the results are stored.
* `configfile`: The path to the config file that contains e.g. the paths to the genome indices (see [Config file](#config-file)).
* `software-deployment-method`: Currently conda and apptainer are supported and tested.
* `workflow-profile`: To define the resources for each snakemake rule.
* `--default-resources`: The temporary directory has to be use by some rules. We advise to set the tmpdir to a fast storage, if available. Otherwise use "tmpdir='<output_directory>'". If not specified, tmpdir is set to the systems default temporary directory.
* `--conda-prefix` (optional): Where should the conda environments be stored. Further information in section [Shared conda prefix](#shared-conda-prefix)
* `--profile` (optional): Path to a profile specification that defines e.g. which executor to use and how many jobs are submitted in parallel (see section [Cluster config](#cluster-execution)).
* `--retries` (optional): Number of retries if a rule fails. When a rule is restarted and the default workflow config is not modified, the RAM for the failed rule is increased, to account for potentially higher RAM requirements.

# Config file

The config file contains variables that configure how the workflow is run. This
includes the paths to the reference and indices.

>Note: Make sure to use absolute paths in the config file.

In the config file the following attributes are specified:

* `sample_sheet`: TSV file with sample paths
* `sra_mode`: Input are SRA identifiers (default: `false`)
* `index_dir`: Directory where the genome lib was build. This directory contains the genome resources and the indices (in future this can be build with [tronmake-genome-lib-builder](https://gitlab.rlp.net/tron/tronmake-genome-lib-builder))
* `bam_input`: Input are (u)BAM files (default: `false`)
* `fraser`: Configuration options for FRASER
    * `min_read`: Minimum number of spliced alignment to consider junction in metric calculation (default: `5`).
    * `mapq_filter`: Minimum MAPQ value to consider read for metric calculation.
* `star`: Configuration options for STAR junctions
    * `min_read`: Minimum number of spliced alignments to consider a junction valid. (default: `5`)
* `requantify`: Configuration options for easyquant and splice2neo
    * `interval_mode`: If set to true, run easyquant in interval mode (default: `true`) 
    * `allow_mismatches`: If set to true, allow mismatches in the junction point area (default: `false`) 
    * `bowtie_k_threshold`: Number of multi-mappers allowed in targeted re-quantification alignment. Setting this to `all` would instruct bowtie to report all possible alignments.
    * `cts_size`: Size of context sequence (+/- bp of exonic sequence). Ideally this should be determined based on the fragment size. (default: `1000`) 

* `chrom-filter`: List of chromosomes for which the splice junctions should be kept. Default are human standard chromosomes.


## Example config file

~~~yaml
{% include "../../../config/config.yaml" %}
~~~