#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Annotate wild type junctions

This script can be used to annotate and select a wild type junction in splice2neo output

Example:
    $ python normalize_star_cpm.py -i STAR/SJ.out.tab -o SJ.normalized_cpm.tsv

@Author: Johannes Hausmann
@Date: 2024-09-25
@Copyright: Copyright 2024, TRON gGmbH, Mainz, Germany
@License: MIT
@Version: 0.0.1
@Status: Development

"""

import pandas as pd
import gffutils
import xxhash
from loguru import logger
from intervaltree import IntervalTree
from collections import defaultdict

ALLOWED_CONTIGS = set(f"chr{i}" for i in range(1,23))
ALLOWED_CONTIGS.add("chrX")
ALLOWED_CONTIGS.add("chrY")

def build_interval_trees(df: pd.DataFrame, db):
    """Build exon interval tree for each transcript

    Args:
        df (_type_): _description_
        db (_type_): _description_

    Returns:
        _type_: _description_
    """
    itree = {}
    transcripts = set(df['tx_id'].values)
    for tx in transcripts:
        tx_exons = [f for f in db.children(tx, featuretype="exon")]
        if tx_exons and not tx in itree:
            for exon in tx_exons:
                itree[tx][exon.start:exon.end+1] = exon
    return itree

def preceed(feature, feature_lst):
    """Get preceeding genomic element of feature

    Given a list of genomic elements, return the preceeding
    element of the given feature. For example: A list of
    exons. [e1, e2, e3, e4, e5].

    Args:
        feature (_type_): _description_
        feature_lst (_type_): _description_

    Returns:
        _type_: _description_
    """
    idx = feature_lst.index(feature)
    # Nothing preceeds the first element of a feature list
    if idx == 0:
        return None, None
    else:
        return idx - 1, feature_lst[idx - 1]

def follows(feature, feature_lst):
    """_summary_

    Args:
        feature (Interval): _description_
        feature_lst (list): _description_

    Returns:
        tuple: Index position and element of following element. Otherwise None
    """
    idx = feature_lst.index(feature)
    # Nothing follows the last element of a list
    if idx == len(feature_lst) - 1:
        return None, None
    else:
        return idx + 1, feature_lst[idx + 1]


def _pos_to_junc(chr: str, start: int, end: int, strand: str) -> list[str]:
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
        if not chr in ALLOWED_CONTIGS:
            raise ValueError("Chromosome annotation not compatible")
        if start > end:
            raise ValueError("Genomic start coordinate can not be greater than end coordinate. Please check the sorting of the GTF file used to align reads")
        
        if strand ==  '*':
            junc_id.extend([f"{chr}:{start}-{end}:+",f"{chr}:{start}-{end}:-"])
        else:
            junc_id.append(f"{chr}:{start}-{end}:{strand}")
        
        return junc_id


def _junc_to_pos(junc_id):
    """junction_id to coordinates

    Args:
        junc_id (str): Splice junction in standardized splice junction format (splice2neo)

    Returns:
        tuple: A tuple with chromosome, start, stop and strand information
    """
    chromosome, coord, strand = junc_id.split(':')
    start, stop = coord.split('-')
    return chromosome, int(start), int(stop), strand
    

def position_in_exon(position, coord):
    """Determine junction position in exon

    Args:
        position (int): Genomic position of splice site
        coord (_type_): Interval object of overlapping exon

    Returns:
        str: A descriptor for the junction position in exon
    """
    start = coord.begin
    stop = coord.end - 1
    if position > start and position < stop:
        return "within"
    elif position == start or position == stop:
        return "boundary"
    else:
        return "outside"

def get_wt_junction(df, itree):
    main_wt_junc = []
    additional_wt_junc = []

    for x in df.itertuples():
        
        tx_id = x.tx_id
        junc_id = x.junc_id
        putative_event_type = x.putative_event_type
        # Get a list of exons from Tree
        tx_exons = sorted(itree[tx_id].items())
        
        chrom, start, stop, strand = _junc_to_pos(junc_id)
        start_canonical = ""
        stop_canonical = ""
        wt_junc_id = None
        wt_juncs = []
        # If there are no exons in transcript --> skip
        if not tx_exons:
            junc_ids_total.append(wt_junc_id)
            continue

        # Find exon overlapping the start
        # If overlapping an exon get boundary position  --> canonical splice site
        overlap_start = sorted(itree[tx_id][start])
        try:
            overlap_start = overlap_start[0]
            boundary_start = position_in_exon(start, overlap_start)
        # Overlaps nothing --> either outside or intron
        except IndexError:
            overlap_start = None
            boundary_start = None
        # Find exon overlapping the stop
        # If overlapping an exon get boundary position --> canonical splice site
        overlap_stop = sorted(itree[tx_id][stop])
        try:
            overlap_stop = overlap_stop[0]
            boundary_stop = position_in_exon(stop, overlap_stop)
        except IndexError:
            overlap_stop = None
            boundary_stop = None
        
        if overlap_start is None and overlap_stop is None:
            logger.warning(f"{junc_id} overlaps no exon features in {tx_id}")
            continue
        
        # Use prior information of splice2neo
        
        if putative_event_type == "ASS":
            # This assumes a canonical splice site on one of the ends
            # A3SS on + and A5SS on -
            if boundary_start == "boundary":
                _, next_exon = follows(overlap_start, tx_exons)
                wt_junc_id = _pos_to_junc(chrom, start, next_exon.begin, strand)
                main_wt_junc.append(wt_junc_id)
                additional_wt_junc.append(None)
            # A5SS on + and A3SS on -
            elif boundary_stop == "boundary":
                _, prev_exon = preceed(overlap_stop, tx_exons)
                wt_junc_id = _pos_to_junc(chrom, prev_exon.end - 1, stop, strand)
                main_wt_junc.append(wt_junc_id)
                additional_wt_junc.append(None)
        # Handle exon skipping events
        elif putative_event_type == "ES":
            # For MultiExonSkip we should do an envelop query and return all exons contained 
            # Return junctions for exon in between the exon skips
            in_between_exons = sorted(itree[tx_id].envelope(start,end))
            # Simple casette exon
            if len(in_between_exons) == 1:
                _, next_exon = follows(overlap_start, tx_exons)
                _, prev_exon = preceed(overlap_stop, tx_exons)
                wt_junc_id_1 = _pos_to_junc(chrom, start, next_exon.begin, strand)
                wt_junc_id_2 = _pos_to_junc(chrom, prev_exon.end - 1, stop, strand)
                main_wt_junc.append(wt_junc_id_1)
                additional_wt_junc.append(wt_junc_id_2)
            # Multiple exon skip
            else:
                multi_skip_junc = []
                exons = [overlap_start].extend(in_between_exons).append(overlap_stop)
                for skipped_exon in range(0, in_between_exons-1):
                    _, next_exon = follows(skipped_exon, tx_exons)
                    multi_skip_junc.append(_pos_to_junc(chrom, skipped_exon.end - 1, next_exon.begin, strand))
                main_wt_junc.append(multi_skip_junc[0])
                additional_wt_junc.extend(multi_skip_junc[1:])

        elif putative_event_type == "exitron":
            # Exitron is located within a single exon
            wt_junc_id = None
        elif putative_event_type == "complex":
            wt_junc_id = None

        junc_ids_total.append(wt_junc_id)
    df['wt_junc_id'] = junc_ids_total
    return df

def process(df, itree):
    df = df.groupby('tx_id', dropna=False).apply(
                lambda sub_df: get_wt_junction(sub_df, itree)
            )
    return df

#input_splice = ""
#logger.info("-> Processing splice junction positions...")
#input_df = pd.read_csv(
#            input_splice,
#            sep="\t",
#            low_memory=False,
#            comment="#")
#input_df = input_df.loc[input_df.putative_event_type != "exitron"]
#db = gffutils.FeatureDB(db_file)

#itree = build_interval_trees(input_df, db)
#process(input_df, itree)