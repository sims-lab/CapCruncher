import glob
import os
import re
import sys
import time
from collections import OrderedDict
from datetime import timedelta
from functools import wraps
from itertools import cycle, groupby
from typing import Callable, IO, Iterable, Union

import click
import pandas as pd
import ujson
import xxhash
from pybedtools import BedTool
from cgatcore.iotools import zap_file
import pybedtools


def open_logfile(fn: str) -> IO:
    '''Handles instances where the log file is sys.stdout'''

    from xopen import xopen

    if not isinstance(fn, type(sys.stdout)):
        return xopen(fn, "w")
    else:
        return fn


def merge_dictionaries(dicts: list) -> dict:
    '''Merges multiple dictionary entries'''
    dict_merged = dict()
    for d in dicts:
        dict_merged.update(d)
    return dict_merged


def invert_dict(d: dict) -> dict:
    '''Inverts key: value pairs into value: key pairs'''
    return {v: k for k, v in d.items()}


def is_on(param: str) -> bool:
    '''Returns True if parameter in "on" values'''
    values = ["true", "t", "on", "yes", "y", "1"]
    if str(param).lower() in values:
        return True


def is_off(param: str):
    '''Returns True if parameter in "off" values'''
    values = ["", "None", "none", "F", "f"]
    if str(param).lower() in values:
        return True


def is_none(param: str) -> bool:
    '''Returns True if parameter is none'''
    values = ["", "none"]
    if str(param).lower() in values:
        return True


def get_human_readable_number_of_bp(bp: int) -> str:
    '''Converts integer into human readable basepair number'''

    if bp < 1000:
        bp = f"{bp}bp"
    elif (bp / 1e3) < 1000:
        bp = f"{bp / 1e3}kb"
    elif (bp / 1e6) < 1000:
        bp = f"{bp / 1e6}mb"

    return bp


def is_valid_bed(bed: Union[str, BedTool]) -> bool:

    '''Returns true if bed file can be opened and has at least 3 columns'''
    try:
        bed = BedTool(bed)
        if bed.field_count(n=1) >= 3:
            return True

    except Exception as e:
        
        if isinstance(e, FileNotFoundError):
            print('Bed file not found')
        
        elif isinstance(e, IndexError):
            print('Wrong number of fields detected, check separator/ number of columns')

        else:
            print(e)
        
        return False
    
    



def bed_has_name(bed: Union[str, BedTool]) -> bool:
    '''Returns true if bed file has at least 4 columns'''
    if isinstance(bed, str):
        bed = BedTool(bed)

    if bed.field_count(n=1) >= 4:
        return True


def bed_has_duplicate_names(bed) -> bool:
    '''Returns true if bed file has no duplicated names'''
    if isinstance(bed, str):
        bed = BedTool(bed)

    df = bed.to_dataframe()
    if not df["name"].duplicated().shape[0] > 1:
        return True


def get_re_site(recognition_site: str = None) -> str:

    """
    Obtains the recogniton sequence for a supplied restriction enzyme or correctly
    formats a supplied recognition sequence.

    Args:
     cut_sequence - DNA sequence to use for fasta digestion e.g. "GATC"
     restriction_enzyme - Name of restriction enzyme e.g. DpnII  (case insensitive)

    Returns:
     recognition sequence e.g. "GATC" 
    
    Raises:
     ValueError: Error if restriction_enzyme is not in known enzymes

    """

    known_enzymes = {
        "dpnii": "GATC",
        "mboi": "GATC",
        "hindiii": "AAGCTT",
        "ecori": "GAATTC",
        "nlaiii": "CATG",
    }

    if re.match(r"[GgAaTtCc]+", recognition_site): # matches a DNA sequence
        cutsite = recognition_site.upper() # Just uppercase convert and return
    
    elif recognition_site.lower() in known_enzymes:
        cutsite = known_enzymes[recognition_site.lower()]
    
    else:
        raise ValueError("No restriction site or recognised enzyme provided")

    return cutsite

def hash_column(col: Iterable, hash_type=64) -> list:
    '''Convinience function to perform hashing using xxhash on an iterable.
       Not vectorised.'''

    hash_dict = {
        32: xxhash.xxh32_intdigest,
        64: xxhash.xxh64_intdigest,
        128: xxhash.xxh128_intdigest,
    }

    hash_func = hash_dict.get(hash_type)

    return [hash_func(v) for v in col]


def split_intervals_on_chrom(intervals: Union[str, BedTool, pd.DataFrame]) -> dict:
    '''Creates dictionary from bed file with the chroms as keys'''

    intervals = convert_bed_to_dataframe(intervals)
    return {chrom: df for chrom, df in intervals.groupby("chrom")}

    


def intersect_bins(bins_1: pd.DataFrame, bins_2: pd.DataFrame, **bedtools_kwargs) -> pd.DataFrame:
    '''Intersects two sets of genomic intervals using bedtools intersect.

       Formats the intersection in a clearer way than pybedtool auto names.
    
    '''

    bt1 = BedTool.from_dataframe(bins_1)
    bt2 = BedTool.from_dataframe(bins_2)
    bt_intersect = bt1.intersect(bt2, **bedtools_kwargs)
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


def load_json(fn) -> dict:
    '''Convinence function to load gziped json file using xopen.'''

    from xopen import xopen

    with xopen(fn) as r:
        d = ujson.load(r)
        return d


def get_timing(task_name=None) -> Callable:
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
    '''Simple class to ensure subcommand order is maintained by click.'''
    def __init__(self, name=None, commands=None, **attrs):
        if commands is None:
            commands = OrderedDict()
        elif not isinstance(commands, OrderedDict):
            commands = OrderedDict(commands)
        click.Group.__init__(self, name=name, commands=commands, **attrs)

    def list_commands(self, ctx):
        return self.commands.keys()

def zap_files(files):
    '''Runs cgatcore zap_files on all inputs'''
    for fn in files:
        zap_file(fn)


def get_ucsc_color(color) -> str:
    '''Converts rgb to UCSC compatable colours'''
    return ",".join([str(int(i * 255)).strip() for i in color])

def get_colors(items: Iterable, colors: Union[Iterable, None] = None):
    '''Extracts the appropriate number of colours for the items
       and formats them for UCSC.'''
    
    import seaborn as sns
    import matplotlib

    if not colors:
        colors = sns.color_palette("rainbow", len(items))
        return [get_ucsc_color(color) for color in colors]
    else:
        colors = [
            matplotlib.colors.to_rgb(color) for color in re.split(r"\s|,|;", colors)
        ]
        return [color for i, color in zip(items, cycle(colors))]

def make_group_track(
    bigwigs: list, key: Union[callable, str, int], overlay=True) -> dict:

    '''Generates a UCSC super track by grouping inputs by the provided key.'''

    import trackhub

    super_tracks_dict = dict()
    for name, bws in groupby(sorted(bigwigs, key=key), key=key):

        bws = list(bws)
        name_sanitized = trackhub.helpers.sanitize(name)
        # Create a super track
        super_track = trackhub.SuperTrack(name=name_sanitized)

        # Create an overlay track
        if overlay:
            overlay_track = trackhub.AggregateTrack(
                name=f"{name_sanitized}_overlay",
                aggregate="transparentOverlay",
                visibility="full",
                tracktype="bigWig",
                maxHeightPixels="8:80:128",
                showSubtrackColorOnUi="on",
                windowingFunction="maximum",
            )

        # Generate entries for all of the tracks for this group
        for bw, color in zip(bws, get_colors(bigwigs)):

            bw_base = (
                os.path.basename(bw)
                .replace(".bigWig", "")
                .replace(".normalised.", "")
            )
            bw_sanitized = trackhub.helpers.sanitize(bw_base)

            track = trackhub.Track(
                name=f"{bw_sanitized}_{name}",
                source=bw,
                visibility="hide",
                color=color,
                autoScale="off",
                tracktype="bigWig",
                windowingFunction="maximum",
            )

            track_sub = trackhub.Track(
                name=f"{bw_sanitized}_{name}_subtrack",
                source=bw,
                visibility="hide",
                color=color,
                autoScale="off",
                tracktype="bigWig",
                windowingFunction="maximum",
            )

            super_track.add_tracks(track)

            if overlay:
                overlay_track.add_subtrack(track_sub)

        if overlay:
            super_track.add_tracks(overlay_track)

        super_tracks_dict[name] = super_track

    return super_tracks_dict


class PysamFakeEntry():
    '''Testing class used to supply a pysam FastqProxy like object'''
    def __init__(self, name, sequence, quality):
        self.name = name
        self.sequence = sequence
        self.quality = quality
        self.comment = ''
    
    def __repr__(self) -> str:
       return  '|'.join([self.name, self.sequence, '+', self.quality])

def convert_to_bedtool(bed: Union[str, BedTool, pd.DataFrame]) -> BedTool:
    '''Converts a str or pd.DataFrame to a pybedtools.BedTool object'''
    if isinstance(bed, str):
        bed_conv = BedTool(bed)
    elif isinstance(bed, pd.DataFrame):
        bed_conv = BedTool.from_dataframe(bed)
    elif isinstance(bed, BedTool):
        bed_conv = bed
    
    return bed_conv

def convert_bed_to_dataframe(bed: Union[str, BedTool, pd.DataFrame]) -> pd.DataFrame:
    '''Converts a bed like object (including paths to bed files) to a pd.DataFrame'''
    
    if isinstance(bed, str):
        bed_conv = BedTool(bed).to_dataframe()
    
    elif isinstance(bed, BedTool):
        bed_conv = bed.to_dataframe()
    
    elif isinstance(bed, pd.DataFrame):
        bed_conv = bed

    
    return bed_conv


def format_coordinates(coordinates: Union[str, os.PathLike]) -> BedTool:
    """Converts coordinates supplied in string format or a .bed file to a BedTool.

    Args:
        coordinates (Union[str, os.PathLike]): Coordinates in the form chr:start-end or a path.
    Raises:
        ValueError: Inputs must be supplied in the correct format.

    Returns:
        BedTool: BedTool object containing the required coordinates.
    """

    coordinates = str(coordinates)
    pattern_genomic_coord = re.compile(r"chr[0-2xXyYmM][0-9]*:\d+-\d+(\s\w)*$")
    pattern_bed_file = re.compile(r"(.*).bed")

    if pattern_genomic_coord.match(coordinates):

        coordinates_split = re.split(":|-", coordinates)
        if len(coordinates_split) < 4:
            coordinates_split.append("region_0")

        bt = BedTool(" ".join(coordinates_split), from_string=True)

    elif pattern_bed_file.match(coordinates) and bed_has_name(coordinates):
        bt = BedTool(coordinates)

    elif pattern_bed_file.match(coordinates) and not bed_has_name(coordinates):
        bt = (
            BedTool(coordinates)
            .to_dataframe()
            .reset_index()
            .assign(name=lambda df: "region_" + df["index"].astype("string"))[
                ["chrom", "start", "end", "name"]
            ]
            .pipe(BedTool.from_dataframe)
        )
    else:
        raise ValueError(
            """Coordinates not provided in the correct format. Provide coordinates in the form chr[NUMBER]:[START]-[END] or a .bed file"""
        )

    return bt

def convert_interval_to_coords(interval: Union[pybedtools.Interval, dict], named=False) -> str:
    """Converts interval object to standard genomic coordinates.

    e.g. chr1:1000-2000 

    Args:
        interval (Union[pybedtools.Interval, dict]): Interval to convert.

    Returns:
        str: Genomic coordinates in the format chr:start-end
    """
    if not named:
        return f'{interval["chrom"]}:{interval["start"]}-{interval["end"]}'
    else:
        return (interval['name'], f'{interval["chrom"]}:{interval["start"]}-{interval["end"]}')