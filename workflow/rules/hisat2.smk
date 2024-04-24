rule hisat2:
    """
    Align reads against reference genome
    """
    input:
        reads = ["results/trimmed/{sample}_R1.fastq.gz", "results/trimmed/{sample}_R2.fastq.gz"],
        idx = get_index,
    output:
        bam = temp("results/hisat2/{sample}/{sample}.bam")
    log:
        "logs/hisat2_align_{sample}.log",
    params:
        extra="-dta",
    threads: 8
    resources:
        mem_mb = 8000 
    wrapper:
        "v3.8.0/bio/hisat2/align"

rule sambamba_sort:
    input:
        rules.hisat2.output.bam
    output:
        bam = "results/hisat2/{sample}/{sample}.sorted.bam"
    params:
        ""  # optional parameters
    log:
        "logs/sambamba-sort/{sample}.log"
    threads: 8
    wrapper:
        "v3.8.0/bio/sambamba/sort"


rule bam2cram:
    input:
        bam = "results/hisat2/{sample}/{sample}.sorted.bam",
        genome = get_genome_fasta
    params:
        out_dir = lambda wildcards, output: os.path.dirname(output.cram)
    output:
        cram = "results/alignment/{sample}/{sample}.cram",
        crai = "results/alignment/{sample}/{sample}.cram.crai",
    conda:
        "../envs/samtools.yaml"
    log: "logs/sam2cram/{sample}.log"
    threads: 1
    shell:
        """
        samtools view -h --reference {input.context_fa} -O CRAM -o {output.cram} {input.bam} &> {log}
        samtools index {output.cram} &>> {log}
        """
