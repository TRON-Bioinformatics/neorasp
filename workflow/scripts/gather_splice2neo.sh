#!/usr/bin/env bash
#
# @Author: Johannes Hausmann
# @Date: 2024-11-01
# @Copyright: Copyright 2024, TRON gGmbH, Mainz, Germany
# @License: MIT
# @Version: 0.0.1
# @Status: Development


#peptide_fasta=(${snakemake_input[peptide_fasta]})
#cat "${peptide_fasta[@]}" > "${snakemake_output[peptide_fasta]}"

IFS=' ' read -r -a peptide_fasta <<< ${snakemake_input[peptide_fasta]}
IFS=' ' read -r -a peptide_junc <<< ${snakemake_input[peptide_junc]}
IFS=' ' read -r -a neofox_annotation <<< ${snakemake_input[neofox_annotation]}

> "${snakemake_output[peptide_fasta]}"
for file in "${peptide_fasta[@]}"; do
    cat "$file" >> "${snakemake_output[peptide_fasta]}"
done

# Combine annotation tables
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