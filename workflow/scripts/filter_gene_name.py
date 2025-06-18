#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Filter junctions

This script can be used to filter splice junction calls
from gene regions matching a user provided regular expressions such
as junctions from HLA, Mitochondria or Immunoglobulin regions.

Example:
    $ python filter_gene_regex.py -i /path/to/input.tsv -o /path/to/output

@Author: Johannes Hausmann
@Date: 2024-09-25
@Copyright: Copyright 2024, TRON gGmbH, Mainz, Germany
@License: MIT
@Version: 0.0.1
@Status: Development

"""
import argparse
import pandas as pd
import numpy as np
import os
from loguru import logger

epilog = "Copyright (c) 2024 TRON gGmbH (See LICENSE for licensing details)"


class FilterGene:
    # Default regex
    exclude_genes_homo_sapiens = {
        "^MT-": "Mitochondrial gene",
        "^HLA-": "HLA gene",
        "^IGH[VDJCG]?": "Immunoglobulin gene",
        "^IGHA[12]": "Immunoglobulin gene",
        "^IGHM": "Immunoglobulin gene",
        "^IGHE": "Immunoglobulin gene",
        "^IGHEP[12]": "Immunoglobulin gene",
        "^IGK[VJC]?": "Immunoglobulin gene",
        "^IGL[VJC]?": "Immunoglobulin gene",
    }
    exclude_genes_mus_musculus = {
        "mt-": "Mitochondrial gene",
        "^H2-": "MHC gene",
        "^Igh[vmdjmgea]?": "Immunoglobulin gene"
    }

    def __init__(self, junction_table, output_path, organism="human", verbose=False) -> None:
        """Parameter initialization"""
        self.df = self._read_tab(junction_table)
        assert all(
            [x in self.df.columns for x in ["junc_id", "hgnc"]]
        ), "File is missing required columns"
        self.output_path = output_path
        self.organism = organism
        self.filter_regex = FilterGene.exclude_genes_homo_sapiens
        if self.organism == "mouse":
            logger.debug(f"-> Organism set to {self.organism}. Applying corresponding gene filters.")
            self.filter_regex = FilterGene.exclude_genes_mus_musculus
            self._log_genes()

    @staticmethod
    def _read_tab(table) -> pd.DataFrame:
        logger.info("-> Reading input table")
        df = pd.read_csv(table, sep="\t")
        return df

    def _log_genes(self) -> None:
        for this_key, this_value in self.filter_regex.items():
            logger.debug(f"-> Applying {this.key} to remove {this.value}")

    def _annotate_with_gene_regex(self):
        """
        Annotate for each Filter regex the match with the
        HGNC column.
        """
        for this_regex, _ in self.filter_regex.items():
            logger.info(f"-> Matching regex: {this_regex} to HGNC column")
            self.df[this_regex] = self.df.hgnc.str.contains(this_regex, regex=True)

    def _summarize_exclude_intention(self):
        logger.info(f"-> Summarizing exlude matches")
        self.df["exclude_gene"] = self.df[list(self.filter_regex.keys())].any(axis=1)

    def run(self):
        self._annotate_with_gene_regex()
        self._summarize_exclude_intention()
        logger.info(f"-> Writing output")
        with open(
            os.path.join(self.output_path, "gene_exclusion_intention.tsv"), "w"
        ) as file_handle:
            intention_df = self.df[["junc_id", "hgnc"] + list(self.filter_regex.keys())]
            intention_df.to_csv(file_handle, index=False, sep="\t")

        with open(
            os.path.join(self.output_path, "sj_problematic_gene.tsv"), "w"
        ) as fail_handle, open(
            os.path.join(self.output_path, "sj_pass_gene.tsv"), "w"
        ) as pass_handle:
            output_df = self.df.drop(list(self.filter_regex.keys()), axis=1)
            output_df[output_df.exclude_gene].to_csv(fail_handle, index=False, sep="\t")
            output_df[~output_df.exclude_gene].to_csv(
                pass_handle, index=False, sep="\t"
            )


def main():
    filter = FilterGene(snakemake.input["parsed_sj"], snakemake.params["working_dir"], organism=snakemake.params["organism"])
    filter.run()


if __name__ == "__main__":
    main()
