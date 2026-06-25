#!/usr/bin/env python3
import os
import pandas as pd


def get_tpm_from_stringie(gffcompare_tmap: str, output_tpm: str):
    """Extract TPM values from StringTie output

    Read gffcompare tmap file and extract TPM values for
    each transfrag assembled by StringTie.
    Args:
        input.gff_tmap (str): Path to gffcompare tmap file
        output.tpm (str): Path to output tsv file with TPM values
    Returns:
        None
    Raises:
        FileNotFoundError: If input file is not found
        ValueError: If input file is not in expected format
    """
    if not os.path.exists(gffcompare_tmap):
        raise FileNotFoundError(f"Input file {gffcompare_tmap} not found")
    df = pd.read_csv(gffcompare_tmap, sep="\t")

    if not all(x in df.columns for x in ["qry_id", "TPM"]):
        raise ValueError("Input file is not in expected format")

    df = df[["qry_id", "TPM"]]
    df.columns = ["stringtie_tx_id", "stringtie_TPM"]

    df.to_csv(output_tpm, sep="\t", index=False)


def junc_to_tpm(junc_to_tx: str, tpm: str, junc_to_tpm: str):
    """Map junctions to TPM values

    Map splice junctions to TPM values based on transcript associations.
    Args:
        junc_to_tx (str): Path to junction to transcript mapping file
        tpm (str): Path to TPM values file
        junc_to_tpm (str): Path to output junction to TPM mapping file
    """
    junc_df = pd.read_csv(junc_to_tx, sep="\t", header=None)
    junc_df.columns = ["chrom", "start", "end", "strand", "tx_ids"]
    tpm_df = pd.read_csv(tpm, sep="\t")
    junc_df["tx_ids"] = junc_df["tx_ids"].str.split(",")
    junc_df = junc_df.explode("tx_ids").reset_index(drop=True)
    merged_df = pd.merge(
        junc_df, tpm_df, left_on="tx_ids", right_on="stringtie_tx_id", how="left"
    )
    merged_df = merged_df[
        ["chrom", "start", "end", "strand", "tx_ids", "stringtie_TPM"]
    ]
    merged_df = (
        merged_df.groupby(["chrom", "start", "end", "strand"])
        .agg({"tx_ids": lambda x: ",".join(x), "stringtie_TPM": "sum"})
        .reset_index()
    )
    merged_df.to_csv(junc_to_tpm, sep="\t", index=False)


if __name__ == "__main__":
    if snakemake.rule == "extract_tpm_from_stringtie":
        get_tpm_from_stringie(snakemake.input.gff_tmap, snakemake.output.tpm)
    elif snakemake.rule == "junc_to_tpm":
        junc_to_tpm(
            snakemake.input.junc_to_tx,
            snakemake.input.tpm,
            snakemake.output.junc_to_tpm,
        )
