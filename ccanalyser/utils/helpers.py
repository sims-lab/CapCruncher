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



def check_files_exist(*args):

    *infiles, outfiles = args


    if isinstance(outfiles, str):
        file_exists = os.path.exists(outfiles)

    elif isinstance(outfiles, (tuple, list)):
        file_exists =  all(os.path.exists(fn) for fn in outfiles)
    
    elif outfiles == None:
        files_exist = False
    

    return (False, 'Output files exists') if file_exists else (True, 'Output files do not exist')


