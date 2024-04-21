#!/usr/bin/env python

def get_final_output():
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
    context_seq = config.get("genome", None)
    return context_seq


def get_genome_annotation(wildcards):
    context_seq = config.get("genome", None)
    return context_seq



def get_aligner_sam(wildcards):
    aligner = config["aligner"]
    if aligner in ["star", "STAR"]:
        return f"results/star/{wildcards.sample}/{wildcards.sample}.cram"
    elif aligner in ["hisat2", "HISAT2"]:
        return f"results/hisat2/{wildcards.sample}/{wildcards.sample}.cram"
    else:
        raise snakemake.WorkflowError("Not a valid aligner selctected. Can not find correct alignment...")
