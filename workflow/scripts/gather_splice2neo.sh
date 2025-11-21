#!/usr/bin/env bash
#
# @Author: Johannes Hausmann
# @Date: 2025-10-14
# @Copyright: Copyright 2025, TRON gGmbH, Mainz, Germany
# @License: MIT
# @Version: 0.0.1
# @Status: Development

exec 2> "${snakemake_log[0]}"

# Combine peptide FASTA files
> "${snakemake_output[peptide_fasta]}"
echo ${snakemake_input[peptide_fasta]}
for file in ${snakemake_input[peptide_fasta]}; do
    cat "${file}" >> "${snakemake_output[peptide_fasta]}"
done

# Combine main splice2neo annotation tables
counter=0
> "${snakemake_output[sj_annot_cts_peptide]}"
for this_file in ${snakemake_input[peptide_junc]}
do

    if (( counter == 0 ))
    then
        cat "${this_file}" > "${snakemake_output[sj_annot_cts_peptide]}"
    else
        tail -n +2 "${this_file}" >> "${snakemake_output[sj_annot_cts_peptide]}"
    fi
    counter=$((counter+1))
done

# Combine neofox annotation tables
counter=0
> "${snakemake_output[neofox_annotation]}"
for this_file in ${snakemake_input[neofox_annotation]}
do

    if (( counter == 0 ))
    then
        cat ${this_file} > "${snakemake_output[neofox_annotation]}"
    else
        tail -n +2 ${this_file} >> "${snakemake_output[neofox_annotation]}"
    fi
    counter=$((counter+1))

done