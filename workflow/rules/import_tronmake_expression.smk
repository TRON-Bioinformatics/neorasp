'''
This script is used to import the TronMake RNA-expression workflow rules.
'''

rule fastp:
    """fastp

    Quality control, adapter trimming and read filtering
    of raw read data using fastp.

    input:
        sample (List[str]): A list with paths to forward and reverse reads.
    output:
        trimmed (List[str]): A list with paths to trimmed forward and reverse reads
        unpaired (str): Path to unpaired/singleton reads.
        html (str): Path to fastp html report.
        json (str): Path to fastp json statistics.
    params:
        extra (str): Additional parameters passed to fastp.

    """
    input:
        sample = lambda wildcards: get_fq(wildcards).values()
    output:
        trimmed = [temp("results/{sample}/fastp/{replicate}/{sample}_R1.fastq.gz"), 
                   temp("results/{sample}/fastp/{replicate}/{sample}_R2.fastq.gz")],
        unpaired = temp("results/{sample}/fastp/{replicate}/{sample}_singletons.fastq.gz"),
        html = "results/{sample}/fastp/{replicate}/{sample}.html",
        json = "results/{sample}/fastp/{replicate}/{sample}.json"
    log:
        "results/{sample}/log/fastp_{replicate}.log"
    params:
        extra = ""
    threads: 2
    container:
        config['container'].get('fastp')
    conda:
        '../envs/fastp.yaml'
    benchmark:
        'results/{sample}/benchmark/fastp_{replicate}_bench.txt'
    shell:
        'fastp '
        '--thread {threads} '
        '--in1 {input.sample[0]} '
        '--in2 {input.sample[1]} '
        '--unpaired1 {output.unpaired} '
        '--unpaired2 {output.unpaired} '
        '--out1 {output.trimmed[0]} '
        '--out2 {output.trimmed[1]} '
        '--html {output.html} '
        '--json {output.json} &> {log}'

rule star:
    """STAR

    Align RNA-seq reads against genome using STAR.
    By default, reads are aligned using the ENCODE3
    pipeline options.

    input
        r1 (str): Path to forward (R1) reads.
        r2 (str): Path to reverse (R2) reads.
    output:
        bam (str): Path to unsorted BAM file.
        log (str): Path to STAR execution log.
        log_final (str): Path to final STAR log file.
        sj (str): Path to high-confidence SJ.out.tab.
        chim_junc (str): Path to chimeric SJ out.tab.
        transcriptome_bam (str): Path to alignment in transcript coordinates.
        forward_wig (str): Coverage of forward strand in bedGraph format.
        reverse_wig (str): Coverage of reverse strand in bedGraph format.
        unmapped_fq1 (str): Path to unmapped forward (R1) reads.
        unmapped_fq2 (str): Path to unmapped reverse (R2) reads.
    params:
        extra (List[str]): Additional parameters passed to STAR.
            Defaults to ENCODE3 options.

    """
    input:
        unpack(get_star_input)
    output:
        bam = temp("results/{sample}/star/Aligned.out.bam"),
        log = "results/{sample}/star/Log.out",
        sj = "results/{sample}/star/SJ.out.tab",
        chim_junc = "results/{sample}/star/Chimeric.out.junction",
        log_final = "results/{sample}/star/Log.final.out",
        transcriptome_bam = temp("results/{sample}/star/Aligned.toTranscriptome.out.bam"),
        unmapped_fq1 = "results/{sample}/star/Unmapped.out.mate1.gz",
        unmapped_fq2 = "results/{sample}/star/Unmapped.out.mate2.gz"
    log:
        "results/{sample}/log/star.log",
    params:
        input_str_fq1 = lambda wildcards, input: ','.join(input.fq1),
        input_str_fq2 = lambda wildcards, input: ','.join(input.fq2),
        # ENCODE3 RNA-seq options
        extra=' '.join(['--outSAMtype BAM Unsorted', 
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
                        '--outReadsUnmapped Fastx']),
        prefix = lambda wildcards, output: os.path.dirname(output.bam),
        index = config['star']['ref'],
        read_cmd =
            lambda wildcards, input: determine_star_read_command(wildcards, input.fq1[0]),
        rg_string = lambda wildcards: get_rg_star, 
    threads: 18
    resources:
        mem_mb = 48000
    container:
        config['container'].get('star')
    conda:
        '../envs/star.yaml'
    benchmark:
        'results/{sample}/benchmark/star_bench.txt'
    shell:
        'STAR '
        '{params.read_cmd} '
        '--runThreadN {threads} '
        '--genomeDir {params.index} '
        '--readFilesIn {params.input_str_fq1} {params.input_str_fq2} '
        '{params.extra} '
        '--outSAMattrRGline {params.rg_string} '
        '--outFileNamePrefix {params.prefix}/ &> {log} ; '
        'test -f {params.prefix}/Unmapped.out.mate1 && gzip {params.prefix}/Unmapped.out.mate1 ; '
        'test -f {params.prefix}/Unmapped.out.mate2 && gzip {params.prefix}/Unmapped.out.mate2'

rule samtools:
    """Samtools

    Create index of BAM file for random access.

    input:
        bam (str): Path to name sorted BAM file.
    output:
        bam (str): 
        bai (str): Path to corresponding BAI index file.

    """
    input:
        bam = "results/{sample}/star/Aligned.out.bam"
    output:
        bam = temp("results/{sample}/star/Aligned.sortedByCoord.out.bam"),
        bai = temp("results/{sample}/star/Aligned.sortedByCoord.out.bam.bai")
    container:
        config['container'].get('samtools')
    conda:
        '../envs/samtools.yaml'
    log:
        "results/{sample}/log/samtools.log"
    threads: 1
    benchmark:
        'results/{sample}/benchmark/samtools_bench.txt'
    shell:
        """
        exec 2> {log}
        samtools sort -o {output.bam} {input.bam}
        samtools index {output.bam}
        """

rule bam2cram:
    """CRAM file

    Create CRAM file of alignment

    input:
        bam (str): Path to coordinate sorted BAM file.
        bai (str): Path to corresponding BAI index file.
    output:
        cram (str): Path to alignment in CRAM format.
        crai (str): Path to corresponding CRAI index file

    """
    input:
        bam = "results/{sample}/star/Aligned.sortedByCoord.out.bam",
        bai = "results/{sample}/star/Aligned.sortedByCoord.out.bam.bai",
        genome = config['reference']['genome']
    output:
        cram = "results/{sample}/star/Aligned.sortedByCoord.out.cram",
        crai = "results/{sample}/star/Aligned.sortedByCoord.out.cram.crai"
    container:
        config['container'].get('samtools')
    conda:
        '../envs/samtools.yaml'
    log:
        'results/{sample}/log/samtools_cram.log'
    threads: 4
    resources:
        mem_mb = 8192
    benchmark:
        'results/{sample}/benchmark/bam2cram_bench.txt'
    shell:
        """
        exec 2> {log}
        samtools view -@ {threads} -m 2G -T {input.genome} -C -o {output.cram} {input.bam}
        samtools index {output.cram}
        """

rule qualimap:
    """QualiMap
    
    Gather quality statistics from RNA-seq read alignment
    using QualiMap. This rule executes the rnaseq subcommand
    of QualiMap.

    input:
        bam (str): Path to coordinate sorted BAM file.
        gtf (str): Path to reference annotation in GTF format.
    output:
        directory (str): Path to qualimap working directoy.
    params:
        java_opts (str): Options for the java virtual machine (JVM).
            Defaults to 'JAVA_OPTS="-Xmx8192M -Djava.awt.headless=true"'
        extra (str): Additional parameters passed to QualiMap.
    
    """
    input:
        bam = rules.samtools.output.bam,
        bai = rules.samtools.output.bai,
        # GTF containing transcript, gene, and exon data
        gtf = config['reference']['annotation']
    output:
        directory("results/{sample}/qualimap")
    log:
        "results/{sample}/log/qualimap.log"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb = 16000
    params:
        java_opts = 'JAVA_OPTS="-Xmx15000M -Djava.awt.headless=true"',
        extra = ""
    container:
        config['container'].get('qualimap')
    conda:
        '../envs/qualimap.yaml'
    benchmark:
        'results/{sample}/benchmark/qualimap_bench.txt'
    shell:
        '{params.java_opts} '
        'qualimap rnaseq {params.extra} '
        '-bam {input.bam} -gtf {input.gtf} '
        '-outdir {output} &> {log}'

rule insert_size:
    """RseQC

    Estimate insert size distribution from aligned RNA-seq reads
    and the reference gene annotation. By default, 10.000.000
    reads-pairs are randomly sampled for the calculation.

    input:
        aln (str): Path to coordinate sorted BAM file.
        refgene (str): Path to reference gene model in BED12 format.
    output:
        reads_inner_distance (str): Path to per-read inner distance table.
        freq (str): Path to per inner distance frequency table.
        pdf (str): Path to pdf graph
        plot_r (str): Path to R script
    params:
        extra (str): Additional parameters passed for inner_distance.py.
            Defaults to '-k 10000000 -q 255'
        out_prefix (str): Path to inner_distance.py working directory.

    """
    input:
        aln = rules.samtools.output.bam,
        bai = rules.samtools.output.bai,
        refgene = config['reference']['annotation_bed']
    output:
        reads_inner_distance = "results/{sample}/metrics/{sample}.inner_distance.txt",
        freq = "results/{sample}/metrics/{sample}.inner_distance_freq.txt",
        pdf = "results/{sample}/metrics/{sample}.inner_distance_plot.pdf",
        plot_r = "results/{sample}/metrics/{sample}.inner_distance_plot.r",
    log:
        'results/{sample}/log/insert_size.log',
    params:
        extra = "-k 10000000 -q 255",
        out_prefix = lambda wildcards, output: output.reads_inner_distance.removesuffix('.inner_distance.txt')
    container:
        config['container'].get('additional_software')
    conda:
        '../envs/rseqc.yaml'
    benchmark:
        'results/{sample}/benchmark/insert_size_bench.txt'
    shell:
        "inner_distance.py " 
        "{params.extra} "
        "--input-file {input.aln} "
        "--refgene {input.refgene} "
        "--out-prefix {params.out_prefix} "
        "&> {log} "

rule junction_saturation:
    """
    RSeQC

    Assess junction saturation by downsampling aligned RNA-seq reads
    and counting detected splice junctions at increasing depths.
    This helps evaluate sequencing depth sufficiency.
    """
    input:
        aln = rules.samtools.output.bam,
        bai = rules.samtools.output.bai,
        refgene = config['reference']['annotation_bed']
    output:
        pdf = "results/{sample}/metrics/{sample}.junctionSaturation_plot.pdf",
        rscript = "results/{sample}/metrics/{sample}.junctionSaturation_plot.r"
    log:
        "results/{sample}/log/junction_saturation.log"
    params:
        extra = "-s 10",  # steps of 10% increments up to 100%
        out_prefix = lambda wildcards, output: output.pdf.removesuffix(".junctionSaturation_plot.pdf")
    container:
        config['container'].get('additional_software')
    conda:
        "../envs/rseqc.yaml"
    resources:
        mem_mb = 20000
    benchmark:
        'results/{sample}/benchmark/junction_saturation_bench.txt'
    shell:
        "junction_saturation.py "
        "{params.extra} "
        "--input-file {input.aln} "
        "--refgene {input.refgene} "
        "--out-prefix {params.out_prefix} "
        "-q 255 -v 5 "
        "&> {log}"

rule read_distribution:
    """
    RSeQC

    Analyze genomic distribution of mapped reads across functional categories
    such as exons, introns, UTRs, promoters, and intergenic regions.
    Requires a BED12-format reference annotation.
    """
    input:
        bam = rules.samtools.output.bam,
        bai = rules.samtools.output.bai,
        refgene = config["reference"]["annotation_bed"]
    output:
        txt = "results/{sample}/metrics/{sample}.read_distribution.txt"
    log:
        "results/{sample}/log/read_distribution.log"
    params:
        extra = ""  # optional additional args (e.g., -l for read length)
    container:
        config['container'].get('additional_software')
    conda:
        "../envs/rseqc.yaml"
    benchmark:
        'results/{sample}/benchmark/read_distribution_bench.txt'
    shell:
        "read_distribution.py "
        "{params.extra} "
        "-i {input.bam} "
        "-r {input.refgene} "
        "> {output.txt} "
        "2> {log}"


rule featurecounts:
    """
    Subread featureCounts

    Quantify aligned RNA-seq reads over annotated gene features.
    Uses exon features grouped by a GTF attribute (e.g. gene_type or gene_id).
    Paired-end mode is enabled by default.
    """
    input:
        bam = rules.samtools.output.bam,
        bai = rules.samtools.output.bai,
        gtf = config['reference']['annotation']
    output:
        counts = "results/{sample}/metrics/{sample}.featureCounts.txt",
        summary = "results/{sample}/metrics/{sample}.featureCounts.txt.summary"
    log:
        "results/{sample}/log/featurecounts.log"
    params:
        threads = 4,
        feature_type = "exon",
        attribute_type = "gene_type",
        prefix = lambda wildcards, output: output.counts.removesuffix(".featureCounts.txt")
    container:
        config['container'].get('featurecounts')
    conda:
        "../envs/subread.yaml"
    threads: 4
    benchmark:
        'results/{sample}/benchmark/featurecounts_bench.txt'
    shell:
        "featureCounts "
        "-T {params.threads} "
        "-a {input.gtf} "
        "-o {params.prefix}.featureCounts.txt "
        "-g {params.attribute_type} "
        "-t {params.feature_type} "
        "-p -B -C "
        "{input.bam} "
        "&> {log}"

#rule tin_score:
#    """
#    RSeQC Transcript integrity estimation
#
#    Estimate Transcript Integrity Number (TIN) scores for each transcript
#    based on RNA-seq coverage across exon regions.
#    """
#    input:
#        bam = rules.samtools.output.bam,
#        bai = rules.samtools.output.bai,
#        gtf = config["reference"]["annotation"]
#    output:
#        tin = "results/{sample}/metrics/{sample}.tin.xls",
#        summary = "results/{sample}/metrics/{sample}.summary.xls"
#    log:
#        "results/{sample}/log/tin_score.log"
#    params:
#        extra = ""
#    container:
#        "docker://tronbioinformatics/tron_data_utils:0.0.1"
#    conda:
#        "../envs/rseqc.yaml"
#    shell:
#        "tin.py "
#        "{params.extra} "
#        "-i {input.bam} "
#        "-r {input.gtf} "
#        "&> {log}"
#
rule salmon:
    """Salmon

    Quantification of gene and transcript expression
    using Salmon. Quantification is performed on STAR
    alignment with transcript coordinates (Aligned.toTranscriptome.out.bam) 
    and reference transcripts (ref_cdna.fa).

    input:
        bam (str): Path to Aligned.toTranscriptome.out.bam
        transcripts (str): Path to reference transcripts in FASTA
    output:
        quant (str): Transcript quantification
        quant_gene (str): Gene quantification (sum of all transcript TPMs)
    params:
        libtype (str): Library type of sequencing reads.
            Defaults to automatic detection = 'A'
        extra (str): Additional parameters passed to salmon execution.
            Default is '--seqBias --gcBias --geneMap'
        outdir (str): Output dirname

    """
    input:
        bam = rules.star.output.transcriptome_bam,
        transcripts = config['reference']['cdna'],
        gtf = config['reference']['annotation']
    params:
        libtype = 'A',
        extra = lambda wildcards, input: f'--seqBias --gcBias --geneMap {input.gtf}',
        outdir = lambda wildcards, output: os.path.dirname(output.quant)
    output:
        quant = 'results/{sample}/salmon_bam/quant.sf',
        quant_gene = 'results/{sample}/salmon_bam/quant.genes.sf'
    resources:
        mem_mb = 16000
    threads: 4
    log:
        'results/{sample}/log/salmon_bam_quant.log'
    container:
        config['container'].get('salmon')
    conda:
        '../envs/salmon.yaml'
    benchmark:
        'results/{sample}/benchmark/salmon_bench.txt'
    shell:
        'salmon quant '
        '--alignments {input.bam} '
        '--targets {input.transcripts} '
        '--libType A '
        '{params.extra} '
        '--threads {threads} '
        '--output {params.outdir} &> {log}'
