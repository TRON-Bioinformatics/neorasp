# Workflow overview

This page provides a short overview of the components used in the pipeline followed by the API of individual rules.

## Workflow components

|    Tool    |                       Link                        |  License  |
| :--------: | :-----------------------------------------------: | :-------: |
|   fastp    |         https://github.com/OpenGene/fastp         |    MIT    |
|    STAR    |         https://github.com/alexdobin/STAR         |    MIT    |
|  samtools  |       https://github.com/samtools/samtools        | MIT/Expat |
|   fraser   |       https://github.com/gagneurlab/FRASER        |    MIT    |
|   salmon   |     https://salmon.readthedocs.io/en/latest/      |  GPL-3.0  |
| easyquant  | https://github.com/TRON-Bioinformatics/easyquant  |    MIT    |
| splice2neo | https://github.com/TRON-Bioinformatics/splice2neo |    MIT    |
| stringtie  |      https://ccb.jhu.edu/software/stringtie/      |    MIT    |

### Detailed pipeline dependencies

SnakeMake comes with integrated package management to retrieve and install all software
required to run the pipeline. The following table gives an overview which conda envrionments or
Docker containers are used by individual steps in the pipeline.

{{ read_table('assets/software.tsv', sep = '\\t') }}

## Workflow rules

This section describes in detail the individual snakemake rules of the workflow,
what input they use, what output they produce and optional parameters.

{% include "assets/docstring.md" %}
