#!/usr/bin/env bash
#
# @Author: Johannes Hausmann
# @Date: 2025-10-14
# @Copyright: Copyright 2025, TRON gGmbH, Mainz, Germany
# @License: MIT
# @Version: 0.0.1
# @Status: Development

# Read input function lists of chunk files into arrays
IFS=' ' read -r -a peptide_fasta <<< ${snakemake_input[peptide_fasta]}
IFS=' ' read -r -a peptide_junc <<< ${snakemake_input[peptide_junc]}
IFS=' ' read -r -a neofox_annotation <<< ${snakemake_input[neofox_annotation]}

# Combine peptide FASTA files
> "${snakemake_output[peptide_fasta]}"
for file in "${peptide_fasta[@]}"; do
    cat "$file" >> "${snakemake_output[peptide_fasta]}"
done

# Combine main splice2neo annotation tables
counter=0
> "${snakemake_output[sj_annot_cts_peptide]}"
for this_file in "${peptide_junc[@]}"
do

    if (( counter == 0 ))
    then
        cat "${this_file}" > "${snakemake_output[sj_annot_cts_peptide]}"
    else
        tail -n +2 "${this_file}" >> "${snakemake_output[sj_annot_cts_peptide]}"
    fi
    ((counter++))

done

# Combine neofox annotation tables
counter=0
> "${snakemake_output[neofox_annotation]}"
for this_file in "${neofox_annotation[@]}"
do

    if (( counter == 0 ))
    then
        cat ${this_file} > "${snakemake_output[neofox_annotation]}"
    else
        tail -n +2 ${this_file} >> "${snakemake_output[neofox_annotation]}"
    fi
    ((counter++))

done