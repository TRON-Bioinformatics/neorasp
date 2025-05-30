#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Normalize splice junction counts

This script can be used normalize the raw splice junction counts produced by
STAR using Counts per Million (CPM).

Example:
    $ python normalize_star_cpm.py -i STAR/SJ.out.tab -o SJ.normalized_cpm.tsv

@Author: Johannes Hausmann
@Date: 2024-09-25
@Copyright: Copyright 2024, TRON gGmbH, Mainz, Germany
@License: MIT
@Version: 0.0.2
@Status: Development

"""
import argparse
import pandas as pd
import numpy as np
from loguru import logger

epilog = "Copyright (c) 2025 TRON gGmbH (See LICENSE for licensing details)"


class CpmNormalization:
    # Standard chromosomes in UCSC notation
    chromosomes = [f"chr{x}" for x in list(range(1, 23)) + ["X", "Y"]]
    # Mapping of STAR codes to intron di-nucleotides
    intron_motifs = {
        0: "non-canonical",
        1: "GT/AG",
        2: "CT/AC",
        3: "GC/AG",
        4: "CT/GC",
        5: "AT/AC",
        6: "GT/AT",
    }
    # Mapping of STAR codes to strand
    strand_encoding = {0: "*", 1: "+", 2: "-"}
    # pandas dtypes for STAR SJ.out.tab
    sj_out_dtypes = {
        "chrom": str,
        "start": np.int32,
        "end": np.int32,
        "strandInfo": np.int32,
        "motifInfo": np.int32,
        "annotInfo": np.int32,
        "uniqueReads": np.int32,
        "mmReads": np.int32,
        "overhang": np.int32,
    }
    # pandas column names for STAR SJ.out.tab
    columns_sj_out = [
        "chrom",
        "start",
        "end",
        "strandInfo",
        "motifInfo",
        "annotInfo",
        "uniqueReads",
        "mmReads",
        "max_overhang",
    ]

    def __init__(self, sj_out_tab, output_file, mapped=True) -> None:
        """Parameter initialization"""
        self.sj_out = sj_out_tab
        self.output_file = output_file
        self.seq_depth = self.get_seq_depth(mapped=mapped)

    def get_seq_depth(self, mapped=False) -> int:
        """Parses a star output log file to get input read counts from the fastq origin"""
        log_file = "{}Log.final.out".format(self.sj_out.rstrip("SJ.out.tab"))
        logger.info("-> Reading STAR Log.final.out to obtain sequencing depth")

        number_of_input_reads = 0
        number_of_uniquely_mapped_reads = 0
        number_of_multimapping_reads = 0
        number_of_chimeric_reads = 0

        with open(log_file, "r") as file_handle:
            for line in file_handle:
                if line.split("|")[0].strip() == "Number of input reads":
                    number_of_input_reads = int(line.split("|")[1].strip())
                elif line.split("|")[0].strip() == "Uniquely mapped reads number":
                    number_of_uniquely_mapped_reads = int(line.split("|")[1].strip())
                elif (
                    line.split("|")[0].strip()
                    == "Number of reads mapped to multiple loci"
                ):
                    number_of_multimapping_reads = int(line.split("|")[1].strip())
                elif line.split("|")[0].strip() == "Number of chimeric reads":
                    number_of_chimeric_reads = int(line.split("|")[1].strip())
                else:
                    continue
        if mapped:
            return (
                number_of_uniquely_mapped_reads
                + number_of_multimapping_reads
                + number_of_chimeric_reads
            )
        else:
            return number_of_input_reads

    def _read_sj_out_tab(self) -> pd.DataFrame:
        """Read STAR SJ.out.tab into DataFrame

        Returns:
            pd.DataFrame: STAR splice junctions
        """
        logger.info("-> Reading STAR SJ.out.tab")
        junctions = pd.read_csv(
            self.sj_out, sep="\t", names=self.columns_sj_out, dtype=self.sj_out_dtypes
        )
        return junctions

    @staticmethod
    def _calc_cpm(read_count: int, seq_depth: int) -> float:
        """CPM calculation

        Calculate count per million (CPM) for a feature (read_count)
        based on sequencing depth.

        Args:
            read_count (int): Number of mapped reads supporting the feature
            seq_depth (int):  Number of total reads in sample.

        Returns:
            float: CPM value
        """
        if seq_depth == 0:
            return read_count
        return (read_count / seq_depth) * 1000000

    def cpm_normalize(self, junction_df: pd.DataFrame, mapped=True) -> pd.DataFrame:
        """CPM normalization

        Normalize counts in junction dataframe with CPM.

        Args:
            junction_df (pd.DataFrame): Splice junctions from STAR SJ.out.tab

        Returns:
            pd.DataFrame: Normalized splice junction counts. Normalization for
                uniquely, multi-mapping and total read counts is appended to DataFrame.
        """
        if self.seq_depth == -1:
            junction_df["jCPM_uniquely_mapped"] = np.nan
            junction_df["jCPM_multi_mapped"] = np.nan
            junction_df["jCPM_total_mapped"] = np.nan
        else:
            junction_df["jCPM_uniquely_mapped"] = junction_df.apply(
                lambda x: self._calc_cpm(x.uniqueReads, self.seq_depth), axis=1
            )
            junction_df["jCPM_multi_mapped"] = junction_df.apply(
                lambda x: self._calc_cpm(x.mmReads, self.seq_depth), axis=1
            )
            junction_df["jCPM_total_mapped"] = junction_df.apply(
                lambda x: self._calc_cpm(x.uniqueReads + x.mmReads, self.seq_depth),
                axis=1,
            )
        return junction_df

    def _pos_to_junc(self, chr: str, start: int, end: int, strand: str) -> list[str]:
        """Generate junc_id

        Generate splice junction id based on genomic coordinates of involved exons.
        If STAR aligner was not able to determine strand based on the intron-motif of
        the spliced alignment, we annotate it as positive and negative strand junction.

        Args:
            chr (str): Chromosome identifier
            start (int): Start position of junction (intron coordinate)
            end (int): End position of junction (intron coordinate)
            strand (str): Strand of junction

        Raises:
            ValueError: If chromosome annotation is not compatible.
            ValueError: If genomic coordinates ascendind sorted

        Returns:
            list[str]: Standardized splice junction id(s).
        """
        junc_id = []
        strand = self.strand_encoding[strand]
        if not chr in self.chromosomes:
            raise ValueError("Chromosome annotation not compatible")
        if start > end:
            raise ValueError(
                "Genomic start coordinate can not be greater than end coordinate. Please check the sorting of the GTF file used to align reads"
            )

        if strand == "*":
            junc_id.append(f"{chr}:{start-1}-{end+1}:+ ; {chr}:{start-1}-{end+1}:-")
        else:
            junc_id.append(f"{chr}:{start-1}-{end+1}:{strand}")

        return junc_id

    def run(self):
        """
        Read SJ.out.tab and Log.final.out in a STAR output directory and normalize splice junctions with CPM.
        """
        junction_df = self._read_sj_out_tab()
        junction_df = junction_df.drop(
            junction_df[~junction_df.chrom.isin(self.chromosomes)].index
        )

        logger.info("-> Generating junc_id")
        junction_df["junc_id"] = junction_df.apply(
            lambda x: self._pos_to_junc(x.chrom, x.start, x.end, x.strandInfo),
            axis=1,
            result_type="expand",
        )
        junction_df = junction_df.explode("junc_id")

        logger.info("-> Obtaining splice site motif")
        junction_df["splice_site_motif"] = junction_df.apply(
            lambda x: self.intron_motifs[x.motifInfo].replace("/", "-"), axis=1
        )

        logger.info("-> Normalizing raw read counts to count per million values")
        junction_df = self.cpm_normalize(junction_df)
        logger.info(f"-> Writing data to file {self.output_file}")

        with open(self.output_file, "w") as file_handle:
            junction_df[
                [
                    "junc_id",
                    "splice_site_motif",
                    "jCPM_uniquely_mapped",
                    "jCPM_multi_mapped",
                    "jCPM_total_mapped",
                ]
            ].to_csv(file_handle, index=False, sep="\t")


def main():

    normalize = CpmNormalization(
        snakemake.input["star_sj"], snakemake.output["star_sj_cpm"], mapped=True
    )
    normalize.run()


if __name__ == "__main__":
    main()
