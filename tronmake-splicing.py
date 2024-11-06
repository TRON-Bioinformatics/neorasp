#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""TronMake RNA-splicing

This script can be used to start the SnakeMake workflow

Example:
    $ python tronmake-splicing.py ....

@Author: Johannes Hausmann
@Date: 2024-11-05
@Copyright: Copyright 2024, TRON gGmbH, Mainz, Germany
@License: MIT
@Version: 0.0.1
@Status: Development

"""

import sys
import argparse
import pathlib
import tempfile
import yaml
import subprocess
from loguru import logger


__version__ = "0.0.1"
__pipeline__ = pathlib.Path(__file__).parent / 'workflow' / 'Snakefile'

epilog = "Copyright (c) 2024 TRON gGmbH (See LICENSE for licensing details)"

def execute_cmd(cmd, working_dir = "."):
    """This function runs a command into a subprocess."""
    logger.info("-> Executing CMD: {}".format(" ".join(cmd)))
    p = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = working_dir, shell=False)
    if p.returncode != 0:
        logger.error(p.stderr)
    return p.returncode


def splicing_pipeline(args):
    wf_config = {}
    wf_config["star"] = {"min_read": args.min_expression}
    wf_config["fraser"] = {"min_read": args.min_expression,
                           "mapq_filter": args.mapq}
    wf_config["samples"] = args.samples
    wf_config["index_dir"] = args.genome_lib
    wf_config["sra_mode"] = False
    wf_config["interleaved_input"] = False
    wf_config["bam_input"] = False
    wf_config["requantify"] = {"interval_mode": True, 
                               "allow_mismatches": False,
                               "bowtie_k_threshold": 200}

    with tempfile.NamedTemporaryFile(mode="w", delete=False, dir=args.workdir) as temp_config:
        yaml.dump(wf_config, temp_config)
        temp_config.close()
        cmd = ['snakemake',
               '--snakefile', str(__pipeline__),
               '--local-cores', str(args.jobs),
               '--jobs', str(args.jobs),
               '--configfile', str(temp_config.name),
               '--use-conda',
               '--directory', str(args.workdir),
               '--rerun-triggers', 'mtime']
        if args.slurm:
            cmd.extend(['--executor', 'slurm'])
        returncode = execute_cmd(cmd)

        if returncode != 0:
            logger.error("-> Command \"{}\" returned non-zero exit status".format(cmd))
            sys.exit(1)
        else:
            logger.info("-> Pipeline finished")

def tronmake_cli():
    parser = argparse.ArgumentParser(
        description="TronMake RNA-splice pipeline v{}".format(__version__),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "--slurm",
        dest="slurm",
        help="Execute snakemake pipeline with slurm support",
        action="store_true"
    )
    parser.add_argument(
        "--samples",
        dest="samples",
        help="Sample sheet (tsv)",
        required=True
    )
    parser.add_argument(
        "--genome_lib",
        dest="genome_lib",
        help="Genome library providing genome annotation and tools indices",
        required=True
    )
    parser.add_argument(
        "--min_expression",
        dest="min_expression",
        help="Minimum number of reads to consider splice junction valid",
        default=5
    )
    parser.add_argument(
        "--mapq",
        dest="mapq",
        help="Use only reads with this MAPQ value for splice junction analysis",
        default=255
    )
    parser.add_argument(
        "--jobs",
        dest="jobs",
        help="Number of local CPUs or number of jobs for slurm submission",
        default=16
    )
    parser.add_argument(
        "--workdir",
        dest="workdir",
        help="Work directory for pipeline execution",
        default=pathlib.Path(__file__).parent
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
