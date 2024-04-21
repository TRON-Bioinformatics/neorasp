rule salmon_quant:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        unpack(get_salmon_reads),
        index=get_salmon_index,
        gtf=get_genome_annotation
    output:
        quant="salmon/{sample}/quant.sf",
        lib="salmon/{sample}/lib_format_counts.json",
    log:
        "logs/salmon/{sample}.log",
    params:
        # optional parameters
        libtype="A",
        extra="",
    threads: 2
    wrapper:
        "v3.8.0/bio/salmon/quant"
