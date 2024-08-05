'''
This script is used to import the fastq2vcf workflow rules.

The pipeline path has to be adapted, otherwise it is not found in the fastq2vcf directory
Unfortunately, if a single param is edited, all other params are set empty, so all parameters
from fastq2vcf have to be mirrored to this file. 
'''


module tronmake_expression:
    # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
    snakefile: os.path.join(workflow.basedir, 'tronmake-rna-expression/workflow/Snakefile')
    # when publicly available, use this instead of submodules
    #snakefile: github("owner/repo", path="workflow/Snakefile", tag="v1.0.0")
    # directly pass the parsed rule to the fastq2vcf workflow
    config: config

use rule * from tronmake_expression as tronmake_expression_*


