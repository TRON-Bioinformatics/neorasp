rule salmon:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz",
        index = get_salmon_index,
        gtf = get_genome_annotation
    output:
        quant = "results/salmon/{sample}/quant.sf",
        lib = "results/salmon/{sample}/lib_format_counts.json",
    log:
        "logs/salmon/{sample}.log",
    params:
        # optional parameters
        libtype="A",
        extra="",
    threads: 2
    wrapper:
        "v3.8.0/bio/salmon/quant"
