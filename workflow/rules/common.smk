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
        final_files.extend(
            collect("results/{sample}/fetchdata/sj_final.tsv", sample = sample.sample_name)
        )
        final_files.extend(
            collect("results/{sample}/fetchdata/sj_final_neofox_annotation.tsv", sample = sample.sample_name)
        )
        final_files.extend(
            collect("results/{sample}/fetchdata/sj_final_peptides.fasta", sample = sample.sample_name)
        )
        final_files.extend(
            collect("results/{sample}/metrics/{sample}.inner_distance.txt", sample = sample.sample_name)
        )
        final_files.extend(
            collect("results/{sample}/metrics/{sample}.junctionSaturation_plot.pdf", sample = sample.sample_name)
        )
        #results/{sample}/metrics/{sample}.summary.xls
        final_files.extend(
           collect("results/{sample}/metrics/{sample}.read_distribution.txt", sample = sample.sample_name)
        )
        final_files.extend(
           collect("results/{sample}/metrics/{sample}.featureCounts.txt", sample = sample.sample_name)
        )
        # Cram files
        final_files.extend(
            collect("results/{sample}/star/Aligned.sortedByCoord.out.cram", sample = sample.sample_name)
        )
    return final_files

def read_sample_sheet(file):
    """
    Read sample sheet with either SRA accessions, fastq or BAM files
    """
    with open(file, "r") as file_handle:
        df = pd.read_csv(file_handle, sep="\t", names=["sample_name", "fq1", "fq2"])
        df["fq1"] = df.fq1.str.split(",")
        df["fq2"] = df.fq2.str.split(",")
        df = df.explode(["fq1","fq2"])
        
        # Check for replicates in input
        # cumcount starts with 0 therefore we add 1
        replicate_df = df.groupby("sample_name")["sample_name"].cumcount()
        replicate_df = replicate_df + 1
        # Merge back to expanded df
        df = pd.concat([df, replicate_df], axis=1).rename(columns={0 : "replicate"})
        df["replicate"] = "Rep" + df["replicate"].astype(str)

    return df

def get_fq(wildcards):
    """
    Get FASTQ files

    Paths are determined based on local or SRA mode
    """
    fq = samples.query('sample_name == @wildcards.sample & replicate == @wildcards.replicate')
    fq1 = fq.get('fq1').item()
    fq2 = fq.get('fq2').item()
    
    return {'r1': fq1, 'r2': fq2}

def get_replicates(wildcards):
    """Return the replicate names for a specific sample.
    """
    # use the sample ID to get the replicate rows in the subsample table and return the replicate names
    return samples.query('sample_name == @wildcards.sample').get('replicate', None).values


def get_star_input(wildcards):
    """Collect the STAR input fastq files.

    STAR takes all replicates of a sample. This function collects all replicates
    for the forward and reverse read and returns them as a dictionary.
    """
    replicates = get_replicates(wildcards)
    return {
        'fq1': [f"results/{wildcards.sample}/fastp/{rep}/{wildcards.sample}_R1.fastq.gz" for rep in replicates],
        'fq2': [f"results/{wildcards.sample}/fastp/{rep}/{wildcards.sample}_R2.fastq.gz" for rep in replicates]
    }

def get_rg_star(wildcards):
    """Get the read group line for STAR, where multiple replicates are handled
    """
    replicates = get_replicates(wildcards)
    return ' , '.join([f'ID:{wildcards.sample}_{rep}\tSM:{wildcards.sample}\tPL:ILLUMINA' for rep in replicates])


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
