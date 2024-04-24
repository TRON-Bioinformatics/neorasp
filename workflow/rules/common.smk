#!/usr/bin/env python

def get_final_output():
    """
    Gather final results
    """
    final_files = []
    for sample in samples.itertuples():
        final_files.extend(
            expand(
                "results/hisat/{sample}/{sample}.cram", sample=sample.sample_name
            )
        )
    return final_files

    
def read_sample_sheet(file):
    """
    Read sample sheet as input
    """
    with open(file, "r") as file_handle:
        samples = pd.read_csv(file_handle,
            sep="\t", names=["sample_name", "fq1", "fq2"])
    return samples


def get_fq(wildcards):
    """
    Get FASTQ files from sample sheets
    """
    fq1 = []
    fq2 = []
    
    fq = samples.query('sample_name == @wildcards.sample')
    fq1 = fq.get('fq1', None).item()
    fq2 = fq.get('fq2', None).item()
    
    return {"fq1": fq1, "fq2": fq2}

def get_index(wildcards):
    context_seq = config.get("hisat", None).get("index")
    return context_seq

def get_genome_fasta(wildcards):
    context_seq = config.get("annotation").get("genome", None)
    return context_seq


def get_genome_annotation(wildcards):
    context_seq = config.get("annotation").get("gtf", None)
    return context_seq
