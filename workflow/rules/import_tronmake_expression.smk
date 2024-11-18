'''
This script is used to import the TronMake RNA-expression workflow rules.
'''

rule get_fastq_SRA:
    """SRA fastq-dump

    Dowload SRA Run accessions from archive. Wildcard
    accession is used by fasterq-dump to retrieve
    FASTQ files

    output:
        fq1 (str): Path to foward reads of accession
        fq2 (str): Path to reverse reads of accession
    """
    output:
        # the wildcard name must be accession, pointing to an SRA number
        fq1 = temp("results/{accession}/sra/{accession}_1.fastq.gz"),
        fq2 = temp("results/{accession}/sra/{accession}_2.fastq.gz")
    params:
        fq1_unzipped = lambda wildcards, output: os.path.splitext(output.fq1)[0],
        fq2_unzipped = lambda wildcards, output: os.path.splitext(output.fq2)[1],
        extra = "--skip-technical --split-3",
        tmpdir = lambda wildcards, output: os.path.dirname(output.fq1),
        mem = ""
    threads: 2
    log:
        'results/{accession}/log/sra_download.log'
    container:
        'docker://quay.io/biocontainers/sra-tools:3.1.1--h4304569_0'
    conda:
        '../envs/ucsc_bedgraph_to_bigwig.yaml'
    shell:
        'fasterq-dump '
        '--temp {params.tmpdir} '
        '--threads {snakemake.threads}'
        '{params.mem} '
        '{params.extra} '
        '{params.tmpdir} ' 
        '{accession} ; '
        'gzip {params.fq1_unzipped} ; '
        'gzip {params.fq2_unzipped} '

rule deinterleave:
    """Deinterleave FASTQ

    When FASTQ input is in interleaved format unwrap the 
    reads into separate temporary forward (R1) and reverse (R2)
    reads.

    input:
        interleaved_fastq (str): Path to interleaved FASTQ file
    output:
        r1 (str): Path to temporary R1 reads.
        r2 (str): Path to temporary R2 reads.

    """
    input:
        interleaved_fastq = get_interleaved_input
    output:
        r1 = temp('results/{sample}/deinterleave/{sample}_R1.fastq'),
        r2 = temp('results/{sample}/deinterleave/{sample}_R2.fastq')
    container:
        'docker://busybox:1.37.0-glibc'
    shell:
        "zcat {input.interleaved_fastq} | "
        "paste - - - - - - - - | "
        "tee >(cut -f 1-4 | tr '\\t' '\\n' > {output.r1}) | "
        "cut -f 5-8 | tr '\\t' '\\n' > {output.r2} "

rule bam2fastq:
    """Convert BAM

    When read input is provided as (uBAM), convert alignment 
    into temporary forward (R1) and reverse (R2) FASTQ files.
    Converting step consists of sorting alignment by query
    name, fixing query name (qname) of alignments (removing _1,_2)
    extension and separation into FASTQ files.

    input:
        bam (str): Path to input BAM file
    output:
        r1 (str): Path to temporary R1 reads.
        r2 (str): Path to temporary R2 reads.       
    params:
        repair_script (str): Path to AWK implementation of qname 
            fixing script.
    
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
        'samtools collate --threads 2 -u -O {input.bam} | awk -f {params.repair_script} | '
        'samtools fastq --threads 2 -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n'

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
        trimmed = ["results/{sample}/fastp/{sample}_R1.fastq.gz", 
                   "results/{sample}/fastp/{sample}_R2.fastq.gz"],
        unpaired = "results/{sample}/fastp/{sample}_singletons.fastq.gz",
        html = "results/{sample}/fastp/{sample}.html",
        json = "results/{sample}/fastp/{sample}.json"
    log:
        "results/{sample}/log/fastp.log"
    params:
        extra = ""
    threads: 2
    container:
        'docker://quay.io/biocontainers/fastp:0.23.4--h125f33a_5'
    conda:
        '../envs/fastp.yaml'
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
        '--json {output.json} '

rule star:
    """STAR

    Align RNA-seq reads against genome using STAR.
    By default, reads are aligned using the ENCODE3
    pipeline options.

    input
        r1 (str): Path to forward (R1) reads.
        r2 (str): Path to reverse (R2) reads.
    output:
        alignment (str): Path to unsorted BAM file.
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
        r1 = "results/{sample}/fastp/{sample}_R1.fastq.gz",
        r2 = "results/{sample}/fastp/{sample}_R2.fastq.gz"
    output:
        alignment = temp("results/{sample}/star/Aligned.out.bam"),
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
    container:
        'docker://quay.io/biocontainers/star:2.7.11b--h43eeafb_2'
    conda:
        '../envs/star.yaml'
    shell:
        'STAR '
        '{params.read_cmd} '
        '--runThreadN {threads} '
        '--genomeDir {params.index} '
        '--readFilesIn {input.r1} {input.r2} '
        '{params.extra} '
        '--outFileNamePrefix {params.prefix}/ &> {log}'

rule bedgraph_to_bigwig:
    """BigWig creation
    
    Convert coverage bedGraph files from STAR to binary
    BigWig for genome wide coverage in IGV.

    input:
        bedGraph_forward (str): Coverage of forward strand in bedGraph format.
        bedGraph_reverse (str): Coverage of reverse strand in bedGraph format.
        chromsizes (str): Path to TSV file describing chromosome sizes.
    output:
        bw_forward (str): Coverage of forward strand in BigWig format.
        bw_reverse (str): Coverage of reverse strand in BigWig format.
    params: 
        extra (str): Additional parameters passed to bedGraphToBigWig.
    """
    input:
        bedGraph_forward = rules.star.output.forward_wig,
        bedGraph_reverse = rules.star.output.reverse_wig,
        chromsizes = os.path.join(config['index_dir'], 'ref_genome.chrom.sizes')
    output:
        bw_forward = "results/{sample}/star/Signal.Unique.str1.bw",
        bw_reverse = "results/{sample}/star/Signal.Unique.str2.bw"
    log:
        "results/{sample}/log/bedgraph2bigwig.log"
    params:
        extra = "" # optional params string
    container:
        'docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:472--h9b8f530_1'
    conda:
        '../envs/ucsc_bedgraph_to_bigwig.yaml'
    shell:
        '''
        bedGraphToBigWig {params.extra} {input.bedGraph_forward} {input.chromsizes} {output.bw_forward} &> {log}
        bedGraphToBigWig {params.extra} {input.bedGraph_reverse} {input.chromsizes} {output.bw_reverse} &>> {log}
        '''

rule samtools:
    """Samtools

    Sort queryname sorted alignment from STAR by coordinate
    and create index for random access.

    input:
        bam (str): Path to queryname sorted BAM file.
    output:
        bam (str): Path to coordinate sorted BAM file.
        bai (str): Path to corresponding BAI index file.

    """
    input:
        bam = "results/{sample}/star/Aligned.out.bam"
    output:
        bam = "results/{sample}/star/Aligned.sortedByCoord.out.bam",
        bai = "results/{sample}/star/Aligned.sortedByCoord.out.bam.bai"
    container:
        'docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0'
    conda:
        '../envs/samtools.yaml'
    threads: 1
    shell:
        """
        samtools sort -@{threads} -m 2G -T ${{tmp}}/SORTtmp_{wildcards.sample} -O bam {input.bam} > {output.bam}
        samtools index {output.bam}
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
        java_opts = 'JAVA_OPTS="-Xmx8192M -Djava.awt.headless=true"',
        extra = ""
    container:
        'docker://quay.io/biocontainers/qualimap:2.3--hdfd78af_0'
    conda:
        '../envs/qualimap.yaml'
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
            Defaults to '-k 10000000'
        out_prefix (str): Path to inner_distance.py working directory.

    """
    input:
        aln = rules.samtools.output.bam,
        refgene = os.path.join(config['index_dir'], 'ref_annot.bed')
    output:
        reads_inner_distance = "results/{sample}/metrics/{sample}.inner_distance.txt",
        freq = "results/{sample}/metrics/{sample}.inner_distance_freq.txt",
        pdf = "results/{sample}/metrics/{sample}.inner_distance_plot.pdf",
        plot_r = "results/{sample}/metrics/{sample}.inner_distance_plot.r",
    log:
        'results/{sample}/log/insert_size.log',
    params:
        extra = "-k 10000000",
        out_prefix = lambda wildcards, ouput: output.reads_inner_distance.removesuffix('.inner_distance.txt')
    container:
        'docker://quay.io/biocontainers/rseqc:5.0.4--pyhdfd78af_0'
    conda:
        '../envs/rseqc.yaml'
    shell:
        "inner_distance.py " 
        "{parmas.extra} "
        "--input-file {input.aln} "
        "--refgene {input.refgene} "
        "--out-prefix {params.out_prefix} "
        "&> {log} "

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
        transcripts = os.path.join(config['index_dir'], 'ref_cdna.fa')
    params:
        libtype = 'A',
        extra = f'--seqBias --gcBias --geneMap {os.path.join(config['index_dir'], 'ref_annot.gtf')}',
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
        'docker://quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2'
    conda:
        '../envs/salmon.yaml'
    shell:
        'salmon quant '
        '--alignments {input.bam} '
        '--targets {input.transcripts} '
        '--libType A '
        '{params.extra} '
        '--threads {threads} '
        '--output {params.outdir} &> {log}'
