'''
This script is used to import the TronMake RNA-expression workflow rules.

The pipeline path has to be adapted, otherwise it is not found in the TronMake directory
Unfortunately, if a single param is edited, all other params are set empty, so all parameters
from TronMake RNA-expression have to be mirrored to this file. 
'''


module tronmake_expression:
    # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
    snakefile: os.path.join(workflow.basedir, 'tronmake-rna-expression/workflow/Snakefile')
    # when publicly available, use this instead of submodules
    #snakefile: github("owner/repo", path="workflow/Snakefile", tag="v1.0.0")
    # directly pass the parsed rule to the fastq2vcf workflow
    config: config

use rule * from tronmake_expression as tronmake_expression_*

use rule salmon_quant_bam from tronmake_expression as tronmake_expression_salmon_quant_bam with:
    params:
        libtype = 'A',
        extra = f'--seqBias --gcBias --geneMap {os.path.join(config['index_dir'], 'ref_annot.gtf')}',
        outdir = lambda wildcards, output: os.path.dirname(output.quant)
    output:
        quant = 'results/{sample}/salmon_bam/quant.sf',
        quant_gene = 'results/{sample}/salmon_bam/quant.genes.sf'
    resources:
        mem_mb = 15000

use rule qualimap from tronmake_expression as tronmake_expression_qualimap with:
    params:
        java_opts="-Xmx8192M -Djava.awt.headless=true"
    resources:
        mem_mb = 10000

        