rule hisat2:
    input:
        reads = get_fq,
        idx = get_index,
    output:
        "results/hisat2/{sample}/{sample}.bam"
    log:
        "logs/hisat2_align_{sample}.log",
    params:
        extra="-dta",
    threads: 8
    resources:
        mem_mb = 8000 
    wrapper:
        "v3.8.0/bio/hisat2/align"

rule bam2cram:
    input:
        bam = get_aligner_sam,
        genome = get_genome_fasta
    params:
        out_dir = lambda wildcards, output: os.path.dirname(output.cram)
    output:
        cram = "results/alignment/{sample}/{sample}.cram",
        crai = "results/alignment/{sample}/{sample}.cram.crai",
    conda:
        "../envs/samtools.yaml"
    log: "logs/sam2cram_{sample}.log"
    threads: 1
    shell:
        """
        samtools view -h --reference {input.context_fa} -O CRAM -o {output.cram} {input.sam} &> {log}
        samtools index {output.cram} &>> {log}
        """
