#!/usr/bin/env python
# coding: utf-8


import argparse
import sys, os
import subprocess
import pandas as pd
import numpy as np
import xxhash
from logzero import logger
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



def main():

    # add options to inputs
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter, description="Add wild-type hits of novel junction sequences.\n"
    )

    parser.add_argument("--input_splice", 
        required=True,
        help="splice2neo annotation table")
    parser.add_argument(
        "--output",
        required=True,
        help="splice2neo table with blast hits of target sequence")
    parser.add_argument("--transcriptome",
        required=True, help="Reference transcriptome to blast against.")
    parser.add_argument("--temp_dir",
        default="/tmp", help="tmp directory")
    parser.add_argument("--mismatch",
        default=0, help="Number of mismatches allowed in blast alignment.")
    parser.add_argument(
        "--threads",
        default=1,
        type=int,
        required=False,
        help="number of threads to use by blast")

    args = parser.parse_args()

    input_splice = args.input_splice
    out_file = args.output
    transcriptome = args.transcriptome
    temp_dir = args.temp_dir
    threads = args.threads


    logger.info("-> Processing splice junction positions...")
    input_df = pd.read_csv(
        input_splice,
        sep="\t",
        low_memory=False,
        comment="#"
    )

    logger.info("-> Extracting context sequence around splice junction...")
    temp = input_df.loc[:, ["cts_id", "cts_seq", "cts_junc_pos"]]
    temp["start"] = temp.apply(lambda x: max(0, x.cts_junc_pos - 20), axis=1)
    temp["stop"] = temp.apply(lambda x: min(x.cts_junc_pos + 20, len(x.cts_seq)), axis=1)
    temp["query_sequence"] = temp.apply(lambda x: x.cts_seq[x.start:x.stop], axis=1)
    temp['query_cts_id'] = temp.apply(lambda x: xxhash.xxh64(x.query_sequence).hexdigest(), axis=1)
    # position - alter_start + neuer
    temp["junc_in_query"] = temp.apply(lambda x: x.cts_junc_pos - (x.start + 1) + 1, axis=1)
    # create the FA file (samtools faidx cant take in stdin)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    faidx_output = os.path.join(temp_dir, "faidx_output.fa")
    sequences = temp[['query_cts_id', 'query_sequence']].drop_duplicates()
    fasta_entries = []
    for row in sequences.itertuples(index=False):
        record = SeqRecord(
            Seq(row.query_sequence),
            id = row.query_cts_id,
            name = row.query_cts_id,
            description='')
        fasta_entries.append(record)
    SeqIO.write(fasta_entries, faidx_output, "fasta")
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Run BLAST against WT transcriptome
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("-> Running blast alignment against transcriptome...")

    # Path to the outfmt6 table
    psl_output = os.path.join(temp_dir, "blast_tx_output.format6")
    # blastn -query test.fa -db ref_genome.fa -max_target_seqs 10000 -outfmt 6 -evalue 1e-10 -num_threads 8 -dust no
    cmd = ["blastn",
           "-query",
           faidx_output,
           "-db",
           transcriptome,
           "-num_threads", 
           str(threads),
           "-max_target_seqs",
           "10000", 
           "-outfmt", 
           "6",
           "-evalue", 
           "1e-10",
           "-dust",
           "no", "-out", psl_output]
    
    if not os.path.exists(os.path.join(temp_dir,"blast_tx_output.ok")):
        subprocess.check_call(cmd)

    subprocess.check_call(["touch", os.path.join(temp_dir,"blast_tx_output.ok")])
    
    # Process the psl output
    logger.info("-> Processing blast output")
    header = [
        "query_cts_id",
        "wt_transcript",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "evalue",
        "bitscore"]

    with open(psl_output, 'r') as file_handle:
        blast_df = pd.read_csv(file_handle, sep='\t', names=header, index_col=False)
    blast_df = blast_df.loc[blast_df.gapopen == 0]
    blast_df = blast_df.loc[blast_df.mismatch <= int(args.mismatch)]
    blast_df = blast_df.loc[blast_df.pident >= 95.0]
    # Merge back to junction dataframe
    temp = temp.join(blast_df.set_index('query_cts_id'), on="query_cts_id", how="left")

    # Check if junction in aligned region --> otherwise skip the alignment
    temp = temp.loc[(temp.junc_in_query > temp.qstart) & (temp.junc_in_query < temp.qend)]
    temp.to_csv(f"{out_file}.tmp", sep="\t", index=False)
    if not len(temp.index) == 0:    
        temp['wt_match'] = temp.apply(lambda x: f'{x.wt_transcript}:{int(x.tstart)}:{int(x.tend)}:{int(x.length)}/{len(x.query_sequence)}', axis = 1)
        temp = temp.loc[:,['cts_id', 'wt_match']].drop_duplicates().groupby('cts_id')['wt_match'].apply(';'.join).reset_index()
        input_df = input_df.join(temp.set_index('cts_id'), on='cts_id')
    else:
        input_df['wt_match'] = np.nan

    input_df.to_csv(out_file, sep="\t", index=False)
    sys.exit(0)
    

if __name__ == "__main__":
    main()