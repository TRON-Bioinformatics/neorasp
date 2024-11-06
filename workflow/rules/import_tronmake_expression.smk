'''
This script is used to import the TronMake RNA-expression workflow rules.
'''

rule get_fastq_SRA:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        fq1 = temp("results/{accession}/sra/{accession}_1.fastq.gz"),
        fq2 = temp("results/{accession}/sra/{accession}_2.fastq.gz")
    params:
        extra="--skip-technical"
    threads: 2
    log:
        'results/{accession}/log/sra_download.log'
    container:
        'docker://quay.io/biocontainers/sra-tools:3.1.1--h4304569_0'
    wrapper:
        "v3.10.2/bio/sra-tools/fasterq-dump"

rule deinterleave:
    """
    When fastq input is interleaved unwrap the fastq files into
    separate temporary forward and reverse reads.
    """
    input:
        interleaved_fastq = get_interleaved_input
    output:
        r1 = temp('results/{sample}/deinterleave/{sample}_R1.fastq'),
        r2 = temp('results/{sample}/deinterleave/{sample}_R2.fastq')
    shell:
        "zcat {input.interleaved_fastq} | "
        "paste - - - - - - - - | "
        "tee >(cut -f 1-4 | tr '\\t' '\\n' > {output.r1}) | "
        "cut -f 5-8 | tr '\\t' '\\n' > {output.r2} "

rule bam2fastq:
    """
    When BAM input is provided, separate reads into fastq files.
    """
    input:
        bam = get_bam_input
    output:
        r1 = temp('results/{sample}/bam2fastq/{sample}_R1.fastq.gz'),
        r2 = temp('results/{sample}/bam2fastq/{sample}_R2.fastq.gz')
    threads: 4
    params: 
        repair_script = workflow.source_path('../scripts/repair_sam_qname.awk')
    container:
        'docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0'
    conda:
        '../envs/samtools.yaml'
    shell:
        'samtools collate --threads 2 -u -O {input.bam} | awk -f {params.repair_script} '
        'samtools fastq --threads 2 -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n'

rule fastp:
    """
    Adapter trimming and read filtering
    """
    input:
        sample = lambda wildcards: get_fq(wildcards).values()
    output:
        trimmed = ["results/{sample}/fastp/{sample}_R1.fastq.gz", 
                   "results/{sample}/fastp/{sample}_R2.fastq.gz"],
        unpaired="results/{sample}/fastp/{sample}_singletons.fastq.gz",
        html="results/{sample}/fastp/{sample}.html",
        json="results/{sample}/fastp/{sample}.json"
    log:
        "results/{sample}/log/fastp.log"
    params:
        extra = ""
    threads: 2
    container:
        'docker://quay.io/biocontainers/fastp:0.23.4--h125f33a_5' 
    wrapper:
        "v3.10.2/bio/fastp"

rule star:
    """
    Align RNA-seq reads using ENCODE3 parameters
    """
    input:
        r1 = "results/{sample}/fastp/{sample}_R1.fastq.gz",
        r2 = "results/{sample}/fastp/{sample}_R2.fastq.gz"
    output:
        alignment = "results/{sample}/star/Aligned.sortedByCoord.out.bam",
        log = "results/{sample}/star/Log.out",
        sj = "results/{sample}/star/SJ.out.tab",
        chim_junc = "results/{sample}/star/Chimeric.out.junction",
        log_final = "results/{sample}/star/Log.final.out",
        transcriptome_bam = temp("results/{sample}/star/Aligned.toTranscriptome.out.bam"),
        forward_wig = "results/{sample}/star/Signal.Unique.str1.out.bg",
        reverse_wig = "results/{sample}/star/Signal.Unique.str2.out.bg",
        unmapped_fq1 = "results/{sample}/star/Unmapped.out.mate1",
        unmapped_fq2 = "results/{sample}/star/Unmapped.out.mate2"
    log:
        "results/{sample}/log/star.log",
    params:
        # ENCODE3 RNA-seq options
        extra=' '.join(['--outSAMtype BAM SortedByCoordinate', 
                        '--outFilterType BySJout',
                        '--alignSJoverhangMin 8',
                        '--alignSJDBoverhangMin 1',
                        '--outFilterMismatchNoverReadLmax 0.04',
                        '--alignIntronMin 20', 
                        '--alignIntronMax 1000000',
                        '--alignMatesGapMax 1000000',
                        '--outSAMstrandField intronMotif',
                        '--chimSegmentMin 20',
                        '--quantMode TranscriptomeSAM',
                        '--outWigType bedGraph ',
                        '--outWigStrand Stranded',
                        '--outReadsUnmapped Fastx',
                        '--limitBAMsortRAM 48000000000']),
        prefix = lambda wildcards, output: os.path.dirname(output.alignment),
        index = os.path.join(config['index_dir'], 'indices', 'star'),
        read_cmd =
            lambda wildcards, input: determine_star_read_command(wildcards, input.r1)
    threads: 18
    resources:
        mem_mb = 48000
    conda:
        '../envs/star.yaml'
    container:
        'docker://quay.io/biocontainers/star:2.7.11b--h43eeafb_2'
    shell:
        'STAR '
        '{params.read_cmd} '
        '--runThreadN {threads} '
        '--genomeDir {params.index} '
        '--readFilesIn {input.r1} {input.r2} '
        '{params.extra} '
        '--outFileNamePrefix {params.prefix}/ &> {log}'

rule bedGraphToBigWig_forward:
    """
    Bigwig with genome wide coverage
    """
    input:
        bedGraph = rules.star.output.forward_wig,
        chromsizes = os.path.join(config['index_dir'], 'ref_genome.chrom.sizes')
    output:
        "results/{sample}/star/Signal.Unique.str1.bw"
    log:
        "results/{sample}/log/bedgraph2bigwig_forward.log"
    params:
        "" # optional params string
    container:
        'docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:455--h2a80c09_1'
    wrapper:
        "v3.11.0/bio/ucsc/bedGraphToBigWig"

rule bedGraphToBigWig_reverse:
    input:
        bedGraph = rules.star.output.reverse_wig,
        chromsizes = os.path.join(config['index_dir'], 'ref_genome.chrom.sizes')
    output:
        "results/{sample}/star/Signal.Unique.str2.bw"
    log:
        "results/{sample}/log/bedgraph2bigwig_reverse.log"
    params:
        "" # optional params string
    container:
        'docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:455--h2a80c09_1'
    wrapper:
        "v3.11.0/bio/ucsc/bedGraphToBigWig"

rule qualimap:
    """
    Gather quality statistics
    """
    input:
        bam = rules.star.output.alignment,
        # GTF containing transcript, gene, and exon data
        gtf = os.path.join(config['index_dir'], 'ref_annot.gtf')
    output:
        directory("results/{sample}/qualimap")
    log:
        "results/{sample}/log/qualimap.log"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb = 8192
    params:
        java_opts = "-Xmx8192M -Djava.awt.headless=true"
    container:
        'docker://quay.io/biocontainers/qualimap:2.3--hdfd78af_0'
    wrapper:
        "v3.10.2/bio/qualimap/rnaseq"

rule insert_size:
    input:
        aln = rules.star.output.transcriptome_bam,
        refgene = os.path.join(config['index_dir'], 'ref_annot.bed')
    output:
        reads_inner_distance = "results/{sample}/metrics/{sample}.inner_distance.txt",
        freq = "results/{sample}/metrics/{sample}.inner_distance_freq.txt",
        pdf = "results/{sample}/metrics/{sample}.inner_distance_plot.pdf",
        plot_r = "results/{sample}/metrics/{sample}.inner_distance_plot.r",
    log:
        'results/{sample}/log/insert_size.log',
    params:
        extra = "-k 10000000"
    wrapper:
        "v4.7.2/bio/rseqc/inner_distance"

rule salmon_quant_bam:
    """
    Rule to quantify RNA expression in TPM using Salmon.
    """
    input:
        bam = rules.star.output.transcriptome_bam,
        transcripts = os.path.join(config['index_dir'], 'ref_cdna.fa')
    params:
        libtype = 'A',
        extra = f'--seqBias --gcBias --geneMap {os.path.join(config['index_dir'], 'ref_annot.gtf')}',
        outdir = lambda wildcards, output: os.path.dirname(output.quant)
    output:
        quant = 'results/{sample}/salmon_bam/quant.sf',
        quant_gene = 'results/{sample}/salmon_bam/quant.genes.sf'
    resources:
        mem_mb = 15000
    threads: 4
    log:
        'results/{sample}/log/salmon_bam_quant.log'
    conda:
        '../envs/salmon.yaml'
    container:
        'docker://quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'
    shell:
        'salmon quant '
        '--alignments {input.bam} '
        '--targets {input.transcripts} '
        '--libType A '
        '{params.extra} '
        '--threads {threads} '
        '--output {params.outdir} &> {log}'
