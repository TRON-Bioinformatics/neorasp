#!/usr/bin/env python

def get_final_output():
    """
    Gather final results
    """
    final_files = []
    for sample in samples.itertuples():
        # Re-quantified novel junctions
        final_files.extend(
            expand(
                "results/{sample}/fetchdata/requantified_sj.tsv", sample=sample.sample_name
            )
        )
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
            collect("results/{sample}/fetchdata/requantified_blast_sj.tsv", sample = sample.sample_name)
        )
    return final_files
