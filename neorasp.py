#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""NeoRasp

This script can be used to start the SnakeMake workflow

Example:
    $ python neorasp.py ....

@Author: Johannes Hausmann
@Date: 2024-11-05
@Copyright: Copyright 2025, TRON gGmbH, Mainz, Germany
@License: MIT
@Version: 0.0.9
@Status: Development

"""

import os
import sys
import argparse
import pathlib
import tempfile
import yaml
import subprocess
from collections import defaultdict
from loguru import logger


__version__ = "0.0.5"
__pipeline__ = pathlib.Path(__file__).parent / "workflow" / "Snakefile"

epilog = "Copyright (c) 2025 TRON gGmbH (See LICENSE for licensing details)"


def convert_paths_to_strings(obj):
    if isinstance(obj, dict):
        return {k: convert_paths_to_strings(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_paths_to_strings(i) for i in obj]
    elif isinstance(obj, pathlib.Path):
        return str(obj)
    else:
        return obj


def get_annotation_files_from_genome_lib(genome_lib):
    """
    Collect all files required to execute the workflow from Tron Genome library folder
    """
    wf_config = defaultdict(dict)
    genome_lib_path = pathlib.Path(genome_lib)
    wf_config["star"] = {"ref": genome_lib_path / "indices" / "star"}
    wf_config["reference"]["genome"] = (
        genome_lib_path / "resources" / "ref_genome.fasta"
    )
    wf_config["reference"]["annotation"] = (
        genome_lib_path / "resources" / "ref_annot.gtf"
    )
    wf_config["reference"]["annotation_bed"] = (
        genome_lib_path / "resources" / "ref_annot.bed"
    )
    wf_config["reference"]["cdna"] = (
        genome_lib_path / "resources" / "ref_transcripts.fasta"
    )
    wf_config["reference"]["chromsizes"] = (
        genome_lib_path / "resources" / "chromosome_sizes.txt"
    )
    wf_config["reference"]["encode_mapability"] = (
        genome_lib_path / "resources" / "mappability" / "encode_exclusion.bed"
    )
    wf_config["reference"]["ucsc_mapability"] = (
        genome_lib_path / "resources" / "mappability" / "ucsc_problematic.bed"
    )
    wf_config["reference"]["ref_transcripts"] = (
        genome_lib_path / "indices" / "R" / "ref_transcripts.Rds"
    )
    wf_config["reference"]["ref_cds"] = (
        genome_lib_path / "indices" / "R" / "ref_cds.Rds"
    )
    wf_config["reference"]["2bit"] = (
        genome_lib_path / "indices" / "R" / "ref_genome.2bit"
    )
    wf_config["reference"]["tx2gene"] = (
        genome_lib_path / "resources" / "ref_annot_transcript2gene.tsv"
    )
    wf_config["reference"]["gene2symbol"] = (
        genome_lib_path / "resources" / "ref_annot_gene2symbol.tsv"
    )
    wf_config["reference"]["canonical_juncs"] = (
        genome_lib_path / "resources" / "ref_annot_splice_sites.tsv"
    )
    wf_config["reference"]["rmsk"] = genome_lib_path / "resources" / "ref_rmsk.Rds"

    return wf_config


def execute_cmd(cmd, working_dir="."):
    """This function runs a command into a subprocess."""
    logger.info("-> Executing CMD: {}".format(" ".join(cmd)))
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=working_dir,
        shell=False,
    )
    if p.returncode != 0:
        logger.error(p.stderr)
    return p.returncode


def find_apptainer_mounts(args) -> set[str]:
    """
    SnakeMake does not support auto-mounting directories from the input sample sheet
    and the genome library into the container. Therefore we do this here by finding
    the parent directories and creating the appropriate apptainer
    """
    genome_lib_path = {
        args.genome_lib,
    }
    input_path_set = set()
    with open(args.samples, "r") as file_handle:
        line = file_handle.read()
        line = line.strip().split("\t")
        fq1, fq2 = line[1:]
        input_path_set.add(os.path.abspath(os.path.dirname(fq1)))
        input_path_set.add(os.path.abspath(os.path.dirname(fq2)))
    return genome_lib_path | input_path_set


def generate_apptainer_mounts(paths: set, mode: str = "ro") -> str:
    """
    Generate apptainer mount commands
    """
    apptainer_mounts = []
    for this_path in paths:
        apptainer_mounts.append(f"--bind {this_path}:{this_path}:{mode}")
    return " ".join(apptainer_mounts)


def splicing_pipeline(args):
    wf_config = dict()
    wf_config["fraser"] = {"min_read": args.min_expression, "mapq_filter": args.mapq}
    wf_config["sample_sheet"] = args.samples
    wf_config["requantify"] = {
        "interval_mode": True,
        "allow_mismatches": False,
        "bowtie_k_threshold": 200,
        "cts_size": 1000,
    }

    wf_config["stringtie"] = {"min_junc_count": 1, "min_junc_anchor": 10}

    wf_config["splice2neo"] = {
        "peptide_flank_size": args.peptide_flank_size,
    }

    wf_config["reliable_calls"] = {
        "min_junction_usage": args.min_psi,
        "min_junction_cpm": args.min_cpm,
    }

    genome_lib_config = get_annotation_files_from_genome_lib(args.genome_lib)
    # wf_config["star"] = wf_config["star"] | dict(genome_lib_config)["star"]
    wf_config["reference"] = dict(genome_lib_config)["reference"]
    input_paths = find_apptainer_mounts(args)
    apptainer_bind_commands = generate_apptainer_mounts(input_paths)

    with tempfile.NamedTemporaryFile(
        mode="w", delete=False, dir=args.workdir
    ) as temp_config:
        yaml.dump(convert_paths_to_strings(dict(wf_config)), temp_config)
        temp_config.close()
        cmd = [
            "snakemake",
            "--snakefile",
            str(__pipeline__),
            "--local-cores",
            str(args.jobs),
            "--jobs",
            str(args.jobs),
            "--configfile",
            str(temp_config.name),
            "--sdm",
            "apptainer",
            "--directory",
            str(args.workdir),
            "--rerun-triggers",
            "mtime",
            "--apptainer-args",
            f"{apptainer_bind_commands}",
        ]
        if args.slurm:
            cmd.extend(["--executor", "slurm"])

        returncode = execute_cmd(cmd)

        if returncode != 0:
            logger.error('-> Command "{}" returned non-zero exit status'.format(cmd))
            sys.exit(1)
        else:
            logger.info("-> Pipeline finished")


def tronmake_cli():
    parser = argparse.ArgumentParser(
        description="NeoRasp pipeline v{}".format(__version__),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "--slurm",
        dest="slurm",
        help="Execute snakemake pipeline with slurm support",
        action="store_true",
    )
    parser.add_argument(
        "--samples", dest="samples", help="Sample sheet (tsv)", required=True
    )
    parser.add_argument(
        "--genome_lib",
        dest="genome_lib",
        help="Genome library providing genome annotation and tools indices",
        required=True,
    )
    parser.add_argument(
        "--min_expression",
        dest="min_expression",
        help="Minimum number of reads to consider splice junction valid",
        default=5,
    )
    parser.add_argument(
        "--mapq",
        dest="mapq",
        help="Use only reads with this MAPQ value for splice junction analysis",
        default=255,
    )
    parser.add_argument(
        "--min_cpm",
        dest="min_cpm",
        help="Junction with minimum of this CPM is considered valid",
        default=0.1,
    )
    parser.add_argument(
        "--min_psi",
        dest="min_psi",
        help="Junction with minimum of this PSI value is considered valid",
        default=0.01,
    )
    parser.add_argument(
        "--jobs",
        dest="jobs",
        help="Number of local CPUs or number of jobs for slurm submission",
        default=16,
    )
    parser.add_argument(
        "--workdir",
        dest="workdir",
        help="Work directory for pipeline execution",
        default=pathlib.Path(__file__).parent,
    )
    parser.add_argument(
        "--peptide_flank_size",
        dest="peptide_flank_size",
        help="Size of flanking sequence around splice junction for peptide generation",
        default=13,
    )

    parser.set_defaults(func=splicing_pipeline)

    args = parser.parse_args()

    try:
        args.func(args)
    except AttributeError as e:
        logger.exception(e)
        parser.parse_args(["--help"])


if __name__ == "__main__":
    tronmake_cli()
