# Usage

> **Information: If possible, use absolute paths when specifying directories or files.**\
> **Recommendation: If possible, configure the workflow with a yaml configfile and not on the command line.**\
> **Apptainer: If possible, use apptainer to run the pipeline. This has been tested. Conda execution is possible but not tested.**

```
snakemake \
    --config sample_sheet=</path/to/samples.tsv> index_dir=</path/to/genome_lib> \
    --directory <output_directory> \
    --software-deployment-method apptainer \
    --configfile <path to my config> \
    --workflow-profile workflow/profiles/default \
    --default-resources "tmpdir='<output_directory>'" \
    [--apptainer-prefix </path/to/apptainer/image/location>] \
    [--apptainer-args '--bind </path/to/mount>'] \
    [--conda-prefix </path/to/conda/env/location>] \
    [--profile </path/to/profile/directory>]
    [--retries <n>]
```

- `config`: After `--config` different parameters can be specified (**NOTE**: the parameters after config are separated by space and do not start with `--`. **Moreover they can be omitted if a configfile is used**):
  - `sample_sheet=`: Define the path to the sample sheet.
  - `index_dir=`: Path to the genome library and genome indices.
- `directory`: Specifies where the results are stored.
- `configfile`: The path to the config file that contains e.g. the paths to the genome indices (see [Config file](#config-file)).
- `software-deployment-method`: Currently conda and apptainer are supported. Apptainer is tested and recommended.
- `workflow-profile` (optional): To define the resources for each snakemake rule. The workflow ships a default profile suitable for most analysis.
- `--default-resources`: The temporary directory has to be use by some rules. We advise to set the tmpdir to a fast storage, if available. Otherwise use "tmpdir='\<output_directory>'". If not specified, tmpdir is set to the systems default temporary directory.
- `--conda-prefix` (optional): Where should the conda environments be stored. Further information in section [Shared conda prefix](#shared-conda-prefix)
- `--apptainer-prefix` (optional): Where should the Apptainer images be stored.
- `--apptainer-args` (required when `--sdm apptainer`): Which directories should be mounted into the container images. This is required when using apptainer.
- `--profile` (optional): Path to a profile specification that defines e.g. which executor to use and how many jobs are submitted in parallel (see section [Cluster config](#cluster-execution)).
- `--retries` (optional): Number of retries if a rule fails. When a rule is restarted and the default workflow config is not modified, the RAM for the failed rule is increased, to account for potentially higher RAM requirements.

## Target Rules

NeoRasp provides multiple target rules to run different parts of the pipeline:

- **Default (`all`)**: Runs the complete pipeline including alignment, junction detection, splice2neo annotation, and MultiQC report.

- **`only_alignment`**: Runs only the alignment step with STAR. Useful for quickly aligning reads and identifying splice junctions without running the full annotation pipeline.

  ```bash
  snakemake --until only_alignment [other options]
  ```

- **`only_splice2neo`**: Runs alignment and splice2neo annotation, but skips quality control metrics and MultiQC report generation. Useful for focusing on junction candidates.

  ```bash
  snakemake --until only_splice2neo [other options]
  ```

## Pathvars for Workflow Reuse

NeoRasp implements Snakemake pathvars to make the workflow more flexible and reusable. Pathvars define generic placeholders for output paths that can be customized when importing the workflow as a module.

Default pathvars:

- `results`: `results/{sample}` - Main results directory
- `logs`: `results/{sample}/log` - Log files
- `benchmarks`: `results/{sample}/benchmark` - Benchmark files
- `multiqc`: `results/report/multiqc` - MultiQC report location

These can be overridden when using NeoRasp as a Snakemake module to integrate with other workflows.

## Config file

The config file contains variables that configure how the workflow is run. This
includes the paths to the reference and indices.

> Note: Make sure to use absolute paths in the config file.

In the config file the following attributes are specified:

- `sample_sheet`: TSV file with sample paths

- `fraser`: Configuration options for FRASER

  - `min_read`: Minimum number of spliced alignment to consider junction in metric calculation (default: `5`).
  - `mapq_filter`: Minimum MAPQ value to consider read for metric calculation.

- `star`: Configuration options for STAR alignment

  - `ref`: Path to STAR reference index directory. Note: The index is now specified as a workflow input rather than a parameter, which improves compatibility with Snakemake storage plugins (e.g., for remote file systems).
  - `extra`: Additional STAR parameters (optional). Defaults to ENCODE3 RNA-seq recommendations.

- `requantify`: Configuration options for easyquant and splice2neo

  - `interval_mode`: If set to true, run easyquant in interval mode (default: `true`)
  - `allow_mismatches`: If set to true, allow mismatches in the junction point area (default: `false`)
  - `bowtie_k_threshold`: Number of multi-mappers allowed in targeted re-quantification alignment. Setting this to `all` would instruct bowtie to report all possible alignments.
  - `cts_size`: Size of context sequence (+/- bp of exonic sequence). Ideally this should be determined based on the fragment size. (default: `1000`)

- `splice2neo`: Configuration options for splice2neo.

  - `peptide_flank_size`: Flanking peptide sequence size (default: `13`). The resulting neoantigen candidate will be $$2\*peptide_flank_size$$
  - `scatter_size`: Size of chunks for scatter-gather execution of splice2neo (default: `1000`). Determines how many junctions are processed in each parallel batch.

- `reliable_calls`: Filter criteria to retain spliced alignments.

  - `min_junction_usage`: Minium Intron Jaccard Index (Splice Usage) (default: `0.05`). Splice junctions with less than 5% usage are discarded
  - `min_junction_cpm`: Minimum CPM normalized uniquely mapped reads to keep junction. (default: `0.1`)

- `stringtie`:

  - `min_junc_count`: Minimum number of reads to keep junction in StringTie assembly (default: `1`)
  - `min_junc_anchor`: Minimum anchor of spliced alignment to keep in assembly (default: `10`)

- `reference`: TronMake Genome Library paths.

  - `genome`: Reference genome in fasta format.
  - `annotation`: Reference annotation in GTF format.
  - `annotation_bed`: Reference annotation in BED12 format.
  - `cdna`: Reference transcript sequences in FASTA format. Must match the reference GTF.
  - `chromsizes`: Chromosome sizes file, e.g. from samtools faidx.
  - `encode_mapability`: ENCODE difficult to map regions.
  - `ucsc_mapability`: UCSC difficult to map regions.
  - `ref_transcripts`: GRangesList of reference transcripts as RDS file.
  - `ref_cds`: GRangesList of reference coding sequence as RDS file.
  - `tx2gene`: A tsv mapping transcripts to genes.
  - `gene2symbol`: A tsv file mapping gene ids to HGNC symbols.
  - `2bit`: Reference genome in 2bit format.
  - `canonical_juncs`: A reference set of canoncial splice junctions.
  - `rmsk`: A GRanges Object of RepeatMasker annotation to identify and flag potenital JETs.

### Example config file

```yaml
{% include "../../../config/config.yaml" %}
```

## Shared conda/apptainer prefix

If you want to collaboratively work with this pipeline, it is helpful to have a shared conda environment and apptainer directory. This allows the specification of `--conda-prefix` / `--apptainer-prefix` to the shared directory. If different users use the pipeline, the same environments are not installed multiple times, saving time and storage. **NOTE: It has to be ensured, that the umask of the pipeline users is u=rwx,g=rwx,o= to allow users of the same group to use the created conda environments properly.**

## Apptainer arguments

Please make sure to mount the appropriate directories into the singularity/apptainer image. By default, snakemake will only mount the current working directory and the location of the workflow into the container. However, if your input (FASTQ files and genome library) is located somewhere else you need to pass the appropriate bind commands to the container. In addition, the `.cache` folder of snakemake must be mounted in the container. Therefore `$HOME/.cache/snakemake` should also be included in the `--bind` command.

`--apptainer-args '--bind /path/to/your/input --bind /home/user/.cache/snakemake'`

## FAQ

**Junction of interest not detected.**

If your junction of interest is not detected it might be filtered out in one of the pipeline steps.

- If the junction was lowly covered it might be filtered out by the reliable call criteria.

- If the junction was thought to be canonical you might find it in: `results/{sample}/fetchdata/detected_sj_canonical.tsv`.

- If the junction falls into hard to map/align regions you might find it in: `results/{sample}/fetchdata/mappability/sj_problematic_mappability.tsv`.

- If your junction falls into a highly polymorphic gene you might find it in: `results/{sample}/fetchdata/splice2neo/gene_name_filter/sj_problematic_gene.tsv`

- If your junction does not overlap a gene you might find it in: `results/{sample}/fetchdata/splice2neo/gene_annot/sj_no_transcript_overlap.tsv`

- If your junctions does not alter the peptide sequence it might be removed by the splice2neo filter `cds_description == "mutated cds"`

**Which HGNC genes are excluded by default**

We exclude by default gene loci from highly polymorphic gene regions such as HLA, T-cell or B-cell receptor genes.
The following regex matches are applied to the HGNC gene ids for filtering.

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
