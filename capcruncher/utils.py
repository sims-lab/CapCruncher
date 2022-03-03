import logging
import os
import pickle
import re
import time
from datetime import timedelta
from functools import wraps
from typing import Callable, Generator, Iterable, List, Literal, Tuple, Union
import numpy as np
import pandas as pd
import pybedtools
import ujson
import xxhash
from pybedtools import BedTool


def read_dataframes(filenames: Iterable, **kwargs):
    dframes = []
    for fn in filenames:
        try:
            df = pd.read_csv(fn, **kwargs)
        except (pd.errors.EmptyDataError):
            logging.warning(f"{fn} is empty")

        if not df.empty:
            dframes.append(df)

    if len(dframes) > 0:
        return dframes
    else:
        raise RuntimeError(
            f"All dataframes supplied are empty or incorrectly formatted: {filenames}"
        )


def invert_dict(d: dict) -> Generator[Tuple[str, str], None, None]:
    """Inverts key: value pairs into value: key pairs"""
    yield from ((v, k) for k, v in d.items())


def is_on(param: str) -> bool:
    """
    Returns True if parameter in "on" values

    On values:
        - true
        - t
        - on
        - yes
        - y
        - 1
    """
    values = ["true", "t", "on", "yes", "y", "1"]
    if str(param).lower() in values:
        return True
    else:
        return False


def is_off(param: str):
    """Returns True if parameter in "off" values"""
    values = ["", "None", "none", "F", "f"]
    if str(param).lower() in values:
        return True
    else:
        return False


def is_none(param: str) -> bool:
    """Returns True if parameter is none"""
    values = ["", "none"]
    if str(param).lower() in values:
        return True
    else:
        return False


def get_human_readable_number_of_bp(bp: int) -> str:
    """Converts integer into human readable basepair number"""

    if bp < 1000:
        bp = f"{bp}bp"
    elif (bp / 1e3) < 1000:
        bp = f"{bp / 1e3}kb"
    elif (bp / 1e6) < 1000:
        bp = f"{bp / 1e6}mb"

    return bp


def is_valid_bed(bed: Union[str, BedTool], verbose=True) -> bool:

    """Returns true if bed file can be opened and has at least 3 columns"""
    try:
        bed = BedTool(bed)
        if bed.field_count(n=1) >= 3:
            return True

    except Exception as e:

        if isinstance(e, FileNotFoundError):
            logging.debug(f"Bed file not found")

        elif isinstance(e, IndexError):
            logging.debug(
                f"Wrong number of fields detected check separator/ number of columns"
            )

        else:
            logging.debug(f"Exception raised {e}")


def bed_has_name(bed: Union[str, BedTool]) -> bool:
    """Returns true if bed file has at least 4 columns"""
    if isinstance(bed, str):
        bed = BedTool(bed)

    if bed.field_count(n=1) >= 4:
        return True


def bed_has_duplicate_names(bed: Union[str, BedTool]) -> bool:
    """Returns true if bed file has no duplicated names"""
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

    if re.match(r"[GgAaTtCc]+", recognition_site):  # matches a DNA sequence
        cutsite = recognition_site.upper()  # Just uppercase convert and return

    elif recognition_site.lower() in known_enzymes:
        cutsite = known_enzymes[recognition_site.lower()]

    else:
        logging.error("No restriction site or recognised enzyme provided")
        raise ValueError("No restriction site or recognised enzyme provided")

    return cutsite


def hash_column(col: Iterable, hash_type=64) -> list:
    """
    Convinience function to perform hashing using xxhash on an iterable.

    Function is **not** vectorised.
    """

    hash_dict = {
        32: xxhash.xxh32_intdigest,
        64: xxhash.xxh64_intdigest,
        128: xxhash.xxh128_intdigest,
    }

    hash_func = hash_dict.get(hash_type)

    return [hash_func(v) for v in col]


def split_intervals_on_chrom(intervals: Union[str, BedTool, pd.DataFrame]) -> dict:
    """Creates dictionary from bed file with the chroms as keys"""

    intervals = convert_bed_to_dataframe(intervals)
    return {chrom: df for chrom, df in intervals.groupby("chrom")}


def intersect_bins(
    bins_1: pd.DataFrame, bins_2: pd.DataFrame, **bedtools_kwargs
) -> pd.DataFrame:
    """Intersects two sets of genomic intervals using bedtools intersect.

    Formats the intersection in a clearer way than pybedtool auto names.

    """

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


def load_dict(fn, format: str, dtype: str = "int") -> dict:
    """Convinence function to load gziped json/pickle file using xopen."""

    import itertools

    from xopen import xopen

    if format == "json":
        with xopen(fn) as r:
            d = ujson.load(r)
    elif format == "pickle":
        with xopen(fn, "rb") as r:
            d = pickle.load(r)

    key_sample = list(itertools.islice(d, 50))
    required_dtype = eval(dtype)

    if all(isinstance(k, required_dtype) for k in key_sample):
        return d
    elif isinstance(d, set):
        return {required_dtype(k) for k in d}
    elif isinstance(d, dict):
        return {
            required_dtype(k): required_dtype(v) if v else None for k, v in d.items()
        }


def save_dict(obj: Union[dict, set], fn: os.PathLike, format: str) -> dict:
    """Convinence function to save [gziped] json/pickle file using xopen."""

    from xopen import xopen

    if format == "json":
        with xopen(fn, "w") as w:
            if isinstance(obj, set):
                d = dict.fromkeys(obj)
            else:
                d = obj
            ujson.dump(d, w)
    elif format == "pickle":
        with xopen(fn, "wb") as w:
            pickle.dump(obj, w)

    return fn


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
            logging.info(f"Completed {task_name} in {time_taken} (hh:mm:ss.ms)")
            return result

        return wrapped

    return wrapper


def convert_to_bedtool(bed: Union[str, BedTool, pd.DataFrame]) -> BedTool:
    """Converts a str or pd.DataFrame to a pybedtools.BedTool object"""
    if isinstance(bed, str):
        bed_conv = BedTool(bed)
    elif isinstance(bed, pd.DataFrame):
        bed_conv = BedTool.from_dataframe(bed)
    elif isinstance(bed, BedTool):
        bed_conv = bed

    return bed_conv


def categorise_tracks(ser: pd.Series) -> list:
    """Gets a series for grouping tracks together

    Args:
        ser (pd.Series): File names to map

    Returns:
        list: Mapping for grouping.
    """
    mapping = {
        "raw": "Replicates",
        "normalised": "Replicates_Scaled",
        "summary": "Samples_Summarised",
        "subtraction": "Samples_Compared",
    }
    categories = []
    for index, value in ser.iteritems():
        for key in mapping:
            if key in value:
                categories.append(mapping[key])

    return categories


def convert_bed_to_dataframe(bed: Union[str, BedTool, pd.DataFrame]) -> pd.DataFrame:
    """Converts a bed like object (including paths to bed files) to a pd.DataFrame"""

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

    elif pattern_bed_file.match(coordinates):
        if is_valid_bed(coordinates):
            if bed_has_name(coordinates):
                bt = BedTool(coordinates)
            else:
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
            raise ValueError("Invalid bed file supplied.")

    else:
        raise ValueError(
            """Coordinates not provided in the correct format. Provide coordinates in the form chr[NUMBER]:[START]-[END] or a .bed file"""
        )

    return bt


def convert_interval_to_coords(
    interval: Union[pybedtools.Interval, dict], named=False
) -> Tuple[str]:
    """Converts interval object to standard genomic coordinates.

    e.g. chr1:1000-2000

    Args:
        interval (Union[pybedtools.Interval, dict]): Interval to convert.

    Returns:
        Tuple: Genomic coordinates in the format chr:start-end
    """
    if not named:
        return (
            "Unnammed",
            f'{interval["chrom"]}:{interval["start"]}-{interval["end"]}',
        )
    else:
        return (
            interval["name"],
            f'{interval["chrom"]}:{interval["start"]}-{interval["end"]}',
        )


def gtf_line_to_bed12_line(df):
    df = df.sort_values(["seqname", "start"])
    geneid = df["geneid"].iloc[0]
    exons = df.query('feature == "exon"')
    chrom = df["seqname"].iloc[0]
    start = str(df["start"].min())
    end = str(df["end"].max())
    strand = df["strand"].iloc[0]
    thick_start = start if strand == "+" else end
    thick_end = thick_start
    color = "0,0,0"
    block_count = str(exons.shape[0])
    block_sizes = ",".join((exons["end"] - exons["start"]).values.astype(str))
    block_starts = ",".join((exons["start"] - int(start)).astype(str))

    return "\t".join(
        [
            chrom,
            start,
            end,
            geneid,
            "0",
            strand,
            thick_start,
            thick_end,
            color,
            block_count,
            block_sizes,
            block_starts,
        ]
    )


def get_file_type(fn: os.PathLike) -> str:
    """
    Determines file type based on extension.

    Args:
        fn (os.PathLike): Path to extract file extension from.

    Returns:
        str: File type
    """

    file_types = {
        "hdf5": "hdf5",
        "hdf": "hdf5",
        "json": "json",
        "tsv": "tsv",
        "h5": "hdf5",
        "pkl": "pickle",
        "pickle": "pickle",
        "parquet": "parquet",
    }

    ext = os.path.splitext(os.path.basename(fn).replace(".gz", ""))[-1].strip(".")

    try:
        return file_types[ext]
    except KeyError as e:
        logging.debug(f"File extension {ext} is not supported")
        raise e


def get_categories_from_hdf5_column(
    path: os.PathLike,
    key: str,
    column: str,
    null_value: Union[Literal["."], int, float] = ".",
) -> List[str]:
    """Extracts all categories from pytables table column.

    Args:
        path (os.PathLike): Path to hdf5 store
        key (str): Table name
        column (str): Column name from which to extract categories
        null_value (Union[Literal[, optional): [description]. Defaults to ".".

    Returns:
        List[str]: Category names
    """

    df_test = pd.read_hdf(path, key, start=0, stop=10)

    # If its a category get from the cat codes
    if isinstance(df_test.dtypes[column], pd.CategoricalDtype):
        return [
            cat
            for cat in df_test[column].cat.categories.values
            if not cat == null_value
        ]

    # Try to extract from the metadata
    try:
        with pd.HDFStore(path, "r") as store:
            # s = store.get_storer(key)
            # values = getattr(s.attrs, column)
            values = store[f"{key}_category_metadata"][column].unique()
            return values
    except AttributeError:
        # Just determine from sampling all of the data (HIGHLY inefficient)
        import dask.dataframe as dd
        import dask.distributed

        with dask.distributed.Client(processes=True) as client:
            values = [
                x
                for x in dd.read_hdf(path, key, columns=column)[column]
                .unique()
                .compute()
            ]
        return values


def get_cooler_uri(store: os.PathLike, viewpoint: str, resolution: Union[str, int]):

    cooler_fragment = r"(?P<store>.*?).hdf5::/(?!.*/resolutions/)(?P<viewpoint>.*?)$"
    cooler_binned = (
        r"(?P<store>.*?).hdf5::/(?P<viewpoint>.*?)/resolutions/(?P<binsize>\d+)$"
    )

    if re.match(cooler_fragment, store):
        if resolution:
            uri = f"{store}/resolutions/{resolution}"
        else:
            uri = store

    elif re.match(cooler_binned, store):
        uri = store

    else:

        if not resolution:
            uri = f"{store}::/{viewpoint}"

        else:
            uri = f"{store}::/{viewpoint}/resolutions/{resolution}"

    return uri


class MockFastqRecord:
    """Testing class used to supply a pysam FastqProxy like object"""

    def __init__(self, name, sequence, quality):
        self.name = name
        self.sequence = sequence
        self.quality = quality
        self.comment = ""

    def __repr__(self) -> str:
        return "|".join([self.name, self.sequence, "+", self.quality])


class MockFastaRecord:
    """Testing class used to supply a pysam FastqProxy like object"""

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def __repr__(self) -> str:
        return f">{self.name}\n{self.sequence}\n"
