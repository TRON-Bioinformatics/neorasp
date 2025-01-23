#!/usr/bin/env python
import magic

def get_final_output():
    """
    Gather final results
    """
    final_files = []
    for sample in samples.itertuples():
        final_files.extend(
            collect("results/{sample}/qualimap", sample = sample.sample_name)
        )
        # BigWig of STAR alignment
        final_files.extend(
            collect("results/{sample}/star/Signal.Unique.str1.bw", sample = sample.sample_name)
        )
        final_files.extend(
            collect("results/{sample}/star/Signal.Unique.str2.bw", sample = sample.sample_name)
        )
        final_files.extend(
            collect("results/{sample}/fetchdata/sj_final.tsv", sample = sample.sample_name)
        )
        final_files.extend(
            collect("results/{sample}/fetchdata/sj_final_neofox_annotation.tsv", sample = sample.sample_name)
        )
        #"results/{sample}/metrics/{sample}.inner_distance.txt",
        final_files.extend(
            collect("results/{sample}/metrics/{sample}.inner_distance.txt", sample = sample.sample_name)
        )
    return final_files

def read_sample_sheet(file):
    """
    Read sample sheet with either SRA accessions, fastq or BAM files
    """
    file_content = []
    sra_mode = config['sra_mode']
    interleaved_fastq = config['interleaved_input']
    bam_input = config['bam_input']
    with open(file, "r") as file_handle:
        if sra_mode:
            samples = pd.read_table(file_handle, names=['sample_name'])
            samples = samples[samples.sample_name.str.startswith('SRR') | samples.sample_name.str.startswith('ERR') | samples.sample_name.str.startswith('DRR')]
        else:
            # When input are interleaved fastq files read table with single fastq
            if interleaved_fastq:
                samples = pd.read_csv(file_handle, sep="\t", names=["sample_name", "fq"])
            elif bam_input:
                samples = pd.read_csv(file_handle, sep="\t", names=["sample_name", "bam"])
            else:
                samples = pd.read_csv(file_handle, sep="\t", names=["sample_name", "fq1", "fq2"])
    return samples

def get_fq(wildcards):
    """
    Get FASTQ files

    Paths are determined based on local or SRA mode
    """
    sample = wildcards.sample
    sra_mode = config['sra_mode']
    # If pipeline is downloading from SRA the paths to the FASTQ files
    # are in the fastp subfolder
    if sra_mode:
        return {'r1': "results/{sample}/sra/{sample}_1.fastq.gz",
                'r2': "results/{sample}/sra/{sample}_2.fastq.gz"}
    
    interleaved_fastq = config['interleaved_input']
    if interleaved_fastq:
        return {'r1': 'results/{sample}/deinterleave/{sample}_R1.fastq',
                'r2': 'results/{sample}/deinterleave/{sample}_R2.fastq'}
    
    fq = samples.query('sample_name == @wildcards.sample')
    fq1 = fq.get('fq1').item()
    fq2 = fq.get('fq2').item()
    
    return {'r1': fq1, 'r2': fq2}

def get_interleaved_input(wildcards):
    """
    Get path to interleaved fastq file
    """
    sample = wildcards.sample
    fq = samples.query('sample_name == @wildcards.sample')
    interleaved_fq = fq.get('fq').item()
    return interleaved_fq

def get_bam_input(wildcards):
    """
    Get path to interleaved fastq file
    """
    sample = wildcards.sample
    bam = samples.query('sample_name == @wildcards.sample')
    bam = fq.get('bam').item()
    return bam

def determine_star_read_command(wildcards, read):
    """
    Determine appropriate read command for STAR alignment.
    To ensure we are using the correct uncompression tool we read the magic byte of the file
    """
    # Default is empty string --> uncompressed FASTQ files

    read_command = ""
    magic_byte = magic.from_file(read)
    if "gzip compressed data" in magic_byte:
        read_command = '--readFilesCommand zcat '
    elif "bzip2 compressed data" in magic_byte:
        read_command = '--readFilesCommand bzcat '
    return read_command