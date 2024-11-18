# Workflow overview

This page provides a quick overview of the components used in the pipeline followed by the API
of individual rules.

## Components of workflow
|Tool | Link | License|
|:----:|:----:|:------:|
|fastp| https://github.com/OpenGene/fastp | MIT|
|STAR | https://github.com/alexdobin/STAR | MIT|
|samtools |https://github.com/samtools/samtools | MIT/Expat|
|fraser	| https://github.com/gagneurlab/FRASER | MIT|
|salmon	| https://salmon.readthedocs.io/en/latest/ | GPL-3.0|
|easyquant |https://github.com/TRON-Bioinformatics/easyquant | MIT|
|splice2neo | https://github.com/TRON-Bioinformatics/splice2neo | MIT|


# Workflow rule

This section describes in detail the individual snakemake rules of the workflow, 
what input they use, what output they produce and optional parameters.

{% include "assets/docstring.md" %}
