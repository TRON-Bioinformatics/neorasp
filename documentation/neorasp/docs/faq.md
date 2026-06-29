## What type of splicing events can be detected?

NeoRasp is focused on the detection of intra-gene exon-exon junctions. While other type of junctions (inter-gene exon-exon junctions and intra-chromosomal junctions) might be found in the intial alignment step they are not considered for further annotation and re-quantification. To detect such type of splicing abberations you should use specialized pipelines such as EasyFuse or CTAT-splicing.

## Can NeoRasp also be used for pre-clinical mouse models

We implemented inital support for mouse RNA-seq analysis in release [*v0.0.6*](https://github.com/TRON-Bioinformatics/neorasp/releases/tag/v0.0.6). However this has not been properly tested and will require more work to run as seamlessly as human analysis.

## Why was my junction of interest not detected?

If your junction of interest is not detected it might be filtered out in one of the pipeline steps.

- If the junction was lowly covered it might be filtered out by the reliable call criteria.

- If the junction was thought to be canonical you might find it in: `results/{sample}/fetchdata/detected_sj_canonical.tsv`.

- If the junction falls into hard to map/align regions you might find it in: `results/{sample}/fetchdata/mappability/sj_problematic_mappability.tsv`.

- If your junction falls into a highly polymorphic gene you might find it in: `results/{sample}/fetchdata/splice2neo/gene_name_filter/sj_problematic_gene.tsv`

- If your junction does not overlap a gene you might find it in: `results/{sample}/fetchdata/splice2neo/gene_annot/sj_no_transcript_overlap.tsv`

- If your junctions does not alter the peptide sequence it might be removed by the splice2neo filter `cds_description == "mutated cds"`

## Which HGNC genes are excluded by default?

We exclude by default gene loci from highly polymorphic gene regions such as HLA, T-cell or B-cell receptor genes.
The following regex matches are applied to the HGNC gene ids for filtering when running the pipeline on human data.

| exclude_gene_pattern | exclude_pattern_intention |
| :------------------: | :-----------------------: |
|         ^MT-         |    Mitochondrial gene     |
|        ^HLA-         |         HLA gene          |
|     ^IGH[VDJCG]?     |    Immunoglobulin gene    |
|      ^IGHA[12]       |    Immunoglobulin gene    |
|        ^IGHM         |    Immunoglobulin gene    |
|        ^IGHE         |    Immunoglobulin gene    |
|      ^IGHEP[12]      |    Immunoglobulin gene    |
|      ^IGK[VJC]?      |    Immunoglobulin gene    |
|      ^IGL[VJC]?      |    Immunoglobulin gene    |
|         ^TRB         |   T cell receptor beta    |

For mouse the following filters are applied:

| exclude_gene_pattern | exclude_pattern_intention |
| :------------------: | :-----------------------: |
|         ^mt-         |    Mitochondrial gene     |
|         ^H2-         |         MHC gene          |
|   ^Igh[vmdjmgea]?    |    Immunoglobulin gene    |

## Why are temporary copies of the R-based BSgenome and TxDb objects created for the splice2neo scatter-gather approach?

These objects are not process-safe. When multiple processes access the same objects on disk simultaneously, we observed incorrect results in transcript annotation and context sequence generation. Creating a temporary copy per process ensures each operation is atomic at the filesystem level, guaranteeing correct results across all processes.

The downside is an increased disk footprint. The splice2neo scatter size (configurable in the main config) is therefore a trade-off between speed and disk usage. When processing large cohorts, the scatter size should not be set too small, as this can generate hundreds or even thousands of temporary copies.
