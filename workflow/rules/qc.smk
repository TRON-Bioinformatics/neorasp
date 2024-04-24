rule fastp:
    """
    Quality control using fastp
    """
    input:
        sample=get_fq
    output:
        trimmed=["results/trimmed/{sample}_R1.fastq.gz",
                          "results/trimmed/{sample}_R2.fastq.gz"],
        unpaired="results/trimmed/{sample}.singletons.fastq.gz",
        html="report/trimmed/{sample}.html",
        json="report/trimmed/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    threads: 2
    wrapper:
        "v3.8.0/bio/fastp"

#rule multiqc:
#    input:
#        "results/"
#    output:
#        "qc/multiqc.{sample}.html",
#        "qc_data/multiqc.{sample}_data.zip",
#    params:
#        extra="--verbose",  # Optional: extra parameters for multiqc.
#    log:
#        "logs/multiqc_{sample}.log",
#    wrapper:
#        "v3.4.1/bio/multiqc"

rule feature_counts:
    input:
        # list of sam or bam files
        samples="results/hisat2/{sample}/{sample}.sorted.bam",
        annotation=get_genome_annotation,
        # optional input
        fasta=get_genome_fasta    # implicitly sets the -G flag
    output:
        multiext(
            "results/featurecounts/{sample}/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
    threads: 2
    params:
        strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra="-O --fracOverlap 0.2 -J -p",
    log:
        "logs/featurecounts/{sample}.log",
    wrapper:
        "v3.8.0/bio/subread/featurecounts"

rule qualimap:
    input:
        # BAM aligned, splicing-aware, to reference genome
        bam="results/hisat2/{sample}/{sample}sorted.bam",
        # GTF containing transcript, gene, and exon data
        gtf=get_genome_annotation
    output:
        directory("results/qualimap/{sample}")
    log:
        "logs/qualimap/{sample}.log"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=4096,
    wrapper:
        "v3.8.0/bio/qualimap/rnaseq"