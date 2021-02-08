import os
import re
from typing import Union

import pandas as pd
import pybedtools
import xxhash
from pybedtools import BedTool
from itertools import combinations
import glob


def is_on(param):
    values = ["true", "t", "on", "yes", "y", "1"]
    if str(param).lower() in values:
        return True


def is_off(param):
    values = ["", "None", "none", "F", "f"]
    if str(param).lower() in values:
        return True


def is_none(param):
    values = ["", "none"]
    if str(param).lower() in values:
        return True


def get_human_readable_number_of_bp(bp: int) -> pd.DataFrame:

    if bp < 1000:
        bp = f"{bp}bp"
    elif (bp / 1e3) < 1000:
        bp = f"{bp / 1e3}kb"
    elif (bp / 1e6) < 1000:
        bp = f"{bp / 1e6}mb"

    return bp


def is_valid_bed(bed):
    try:
        bed = BedTool(bed)
        if bed.field_count(n=1) >= 3:
            return True
    except FileNotFoundError:
        return False


def bed_has_name(bed):
    if isinstance(bed, str):
        bed = BedTool(bed)

    if bed.field_count(n=1) >= 4:
        return True


def bed_has_duplicate_names(bed):
    if isinstance(bed, str):
        bed = BedTool(bed)

    df = bed.to_dataframe()
    if not df["name"].duplicated().shape[0] > 1:
        return True


def get_re_site(recognition_site=None):

    """
    Obtains the recogniton sequence for a supplied restriction enzyme or correctly
    formats a supplied recognition sequence.

    Args:
        cut_sequence - DNA sequence to use for fasta digestion e.g. "GATC"
        restriction_enzyme - Name of restriction enzyme e.g. DpnII  (case insensitive)

    Returns:
        recognition sequence e.g. "GATC"

    Raises:
        ValueError if restriction_enzyme is not in known enzymes

    """

    known_enzymes = {
        "dpnii": "GATC",
        "mboi": "GATC",
        "hindiii": "AAGCTT",
        "ecori": "GAATTC",
        "nlaiii": "CATG",
    }

    if re.match(r"[GgAaTtCc]+", recognition_site):
        # This matches a DNA sequence so just convert to upper case and return
        return recognition_site.upper()
    elif recognition_site.lower() in known_enzymes:
        return known_enzymes.get(recognition_site.lower())
    else:
        raise ValueError("No restriction site or recognised enzyme provided")


def collate_histogram_data(fnames):
    return (
        pd.concat([pd.read_csv(fn) for fn in fnames])
        .groupby(["sample", "read_type", "read_number", "number_of_slices"])["count"]
        .sum()
        .reset_index()
        .sort_values(["sample", "read_type", "number_of_slices"])
    )


def collate_read_data(fnames):
    return (
        pd.concat([pd.read_csv(fn) for fn in fnames])
        .groupby(["sample", "stage", "read_type", "read_number", "stat_type"])["stat"]
        .sum()
        .reset_index()
        .sort_values(["sample", "stat"], ascending=[True, False])
    )

def collate_slice_data(fnames):

    df = pd.concat([pd.read_csv(fn) for fn in fnames])
    aggregations = {col: 'sum' if not 'unique' in col else 'max' for col in df.columns
                    if not col in ['sample', 'stage', 'read_type']}

    return (df.groupby(["sample", "stage", "read_type"])
              .agg(aggregations)
              .reset_index()
              .sort_values(["sample", "unique_slices"], ascending=[True, False])
        )

def collate_cis_trans_data(fnames):

    return (pd.concat([pd.read_csv(fn) for fn in fnames])
              .groupby(['sample', 'capture', 'read_type', 'cis/trans'])
              .sum()
              .reset_index()
              .sort_values(["sample", "read_type", 'count'], ascending=[True, True, False]))


def hash_column(col, hash_type=64):

    hash_dict = {
        32: xxhash.xxh32_intdigest,
        64: xxhash.xxh64_intdigest,
        128: xxhash.xxh128_intdigest,
    }

    hash_func = hash_dict.get(hash_type)

    return [hash_func(v) for v in col]

def extract_trimming_stats(fn):
    stat_regexes = {
                'reads_total': re.compile(r'^Total reads processed:\s+([0-9,]+)$'),
                'adapters_removed': re.compile(r'Reads with adapters:\s+([0-9,]+).*'),
                'reads_after_filtering': re.compile(r'Reads written \(passing filters\):\s+([0-9,]+).*'),}
    
    sample_re_match = re.match(r'.*/(.*)_part\d+_(1|2).*', fn)
    
    stats = {}
    stats['sample'] = sample_re_match.group(1)
    stats['read_number'] = sample_re_match.group(2)
    stats['read_type'] = 'pe'
    
    
    with open(fn) as r:
        for line in r:
            for stat_name, pattern in stat_regexes.items():
                regex_match = pattern.match(line)
                if regex_match:
                    stats[stat_name] = int(regex_match.group(1).replace(',', ''))
    
    stats['reads_filtered'] = stats['reads_total'] - stats['reads_after_filtering']
    
    return stats
        
    
def split_intervals_on_chrom(intervals):
    if isinstance(intervals, BedTool):
        intervals = BedTool.to_dataframe()
    elif isinstance(intervals, str):
        intervals = BedTool(intervals).to_dataframe()
        
    return {chrom: df for chrom, df in intervals.groupby('chrom')}

def intersect_bins(bins_1, bins_2):
        bt1 = BedTool.from_dataframe(bins_1)
        bt2 = BedTool.from_dataframe(bins_2)
        bt_intersect = bt1.intersect(bt2, wao=True, sorted=True)
        df_intersect =  (bt_intersect.to_dataframe(disable_auto_names=True, 
                                              header=None,
                                              index_col=False,
                                              names=['chrom_1', 'start_1', 'end_1', 'name_1',
                                                     'chrom_2', 'start_2', 'end_2', 'name_2',
                                                     'overlap'
                                                    ]))
        
        return df_intersect
