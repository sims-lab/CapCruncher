import glob
import os
import re
import sys
import time
from collections import OrderedDict
from datetime import timedelta
from functools import wraps
from itertools import combinations
from typing import Union

import click
import pandas as pd
import pybedtools
import ujson
import xxhash
from pybedtools import BedTool
from xopen import xopen


def open_logfile(fn):
    if not isinstance(fn, type(sys.stdout)):
        return xopen(fn, "w")
    else:
        return fn


def merge_dictionaries(dicts: list):
    dict_merged = dict()
    for d in dicts:
        dict_merged.update(d)
    return dict_merged


def invert_dict(d):
    return {v: k for k, v in d.items()}


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


def hash_column(col, hash_type=64):

    hash_dict = {
        32: xxhash.xxh32_intdigest,
        64: xxhash.xxh64_intdigest,
        128: xxhash.xxh128_intdigest,
    }

    hash_func = hash_dict.get(hash_type)

    return [hash_func(v) for v in col]


def split_intervals_on_chrom(intervals):
    if isinstance(intervals, BedTool):
        intervals = BedTool.to_dataframe()
    elif isinstance(intervals, str):
        intervals = BedTool(intervals).to_dataframe()

    return {chrom: df for chrom, df in intervals.groupby("chrom")}


def intersect_bins(bins_1, bins_2):
    bt1 = BedTool.from_dataframe(bins_1)
    bt2 = BedTool.from_dataframe(bins_2)
    bt_intersect = bt1.intersect(bt2, wao=True, sorted=True)
    df_intersect = bt_intersect.to_dataframe(
        disable_auto_names=True,
        header=None,
        index_col=False,
        names=[
            "chrom_1",
            "start_1",
            "end_1",
            "name_1",
            "chrom_2",
            "start_2",
            "end_2",
            "name_2",
            "overlap",
        ],
    )

    return df_intersect


def load_json(fn):
    with xopen(fn) as r:
        d = ujson.load(r)
        return d


def get_timing(task_name=None):
    """Decorator:
    Gets the time taken by the wrapped function
    """

    def wrapper(f):
        @wraps(f)
        def wrapped(*args, **kwargs):
            time_start = time.perf_counter()
            result = f(*args, **kwargs)
            time_end = time.perf_counter()

            time_taken = timedelta(seconds=(time_end - time_start))
            print(f"Completed {task_name} in {time_taken} (hh:mm:ss.ms)")
            return result

        return wrapped

    return wrapper


class NaturalOrderGroup(click.Group):
    def __init__(self, name=None, commands=None, **attrs):
        if commands is None:
            commands = OrderedDict()
        elif not isinstance(commands, OrderedDict):
            commands = OrderedDict(commands)
        click.Group.__init__(self, name=name, commands=commands, **attrs)

    def list_commands(self, ctx):
        return self.commands.keys()
