import itertools
import os
import pickle
import re
from pathlib import Path
from functools import wraps
from typing import Callable, Iterable, Tuple, Union

import pandas as pd
import pyranges1 as pr


def cycle_argument(arg):
    """Allows for the same argument to be stated once but repeated for all files"""

    if len(arg) == 1:
        return itertools.cycle((arg[0],))
    else:
        return arg


def read_dataframes(filenames: Iterable, **kwargs):
    from loguru import logger

    dframes = []
    for fn in filenames:
        try:
            df = pd.read_csv(fn, **kwargs)
        except pd.errors.EmptyDataError:
            logger.warning(f"{fn} is empty")

        if not df.empty:
            dframes.append(df)

    if len(dframes) > 0:
        return dframes
    else:
        raise RuntimeError(
            f"All dataframes supplied are empty or incorrectly formatted: {filenames}"
        )


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
    return str(param).lower() in values


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


def _read_bed_dataframe(
    bed: Union[str, os.PathLike, pd.DataFrame, pr.PyRanges], nrows=None
) -> pd.DataFrame:
    if isinstance(bed, pd.DataFrame):
        return bed.copy()

    if isinstance(bed, pr.PyRanges):
        return bed.copy()

    df = pd.read_csv(
        bed,
        sep="\t",
        header=None,
        comment="#",
        nrows=nrows,
    )

    column_names = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thick_start",
        "thick_end",
        "item_rgb",
        "block_count",
        "block_sizes",
        "block_starts",
    ]
    df.columns = [
        column_names[i] if i < len(column_names) else f"col_{i}"
        for i in range(df.shape[1])
    ]
    return df


def _standardize_bed_columns(
    df: pd.DataFrame, capitalized: bool = False
) -> pd.DataFrame:
    rename_map = {}
    for column in df.columns:
        column_str = str(column)
        lower = column_str.lower()
        if lower in {"chrom", "chromosome"}:
            rename_map[column] = "Chromosome" if capitalized else "chrom"
        elif lower == "start":
            rename_map[column] = "Start" if capitalized else "start"
        elif lower == "end":
            rename_map[column] = "End" if capitalized else "end"
        elif lower == "name":
            rename_map[column] = "Name" if capitalized else "name"
        elif lower == "score":
            rename_map[column] = "Score" if capitalized else "score"
        elif lower == "strand":
            rename_map[column] = "Strand" if capitalized else "strand"
        elif lower == "thick_start":
            rename_map[column] = "ThickStart" if capitalized else "thick_start"
        elif lower == "thick_end":
            rename_map[column] = "ThickEnd" if capitalized else "thick_end"
        elif lower == "item_rgb":
            rename_map[column] = "ItemRGB" if capitalized else "item_rgb"
        elif lower == "block_count":
            rename_map[column] = "BlockCount" if capitalized else "block_count"
        elif lower == "block_sizes":
            rename_map[column] = "BlockSizes" if capitalized else "block_sizes"
        elif lower == "block_starts":
            rename_map[column] = "BlockStarts" if capitalized else "block_starts"

    return df.rename(columns=rename_map)


def _prepare_intersection_frame(
    df: Union[str, os.PathLike, pd.DataFrame, pr.PyRanges], name_prefix: str
) -> pd.DataFrame:
    frame = convert_bed_to_dataframe(df)
    if frame.empty:
        return frame

    frame = _standardize_bed_columns(frame, capitalized=False)
    for column in ("start", "end"):
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="raise")
    if "name" not in frame.columns:
        frame = frame.copy()
        frame["name"] = [f"{name_prefix}_{idx}" for idx in range(frame.shape[0])]
    else:
        frame["name"] = frame["name"].fillna(
            pd.Series(
                [f"{name_prefix}_{idx}" for idx in range(frame.shape[0])],
                index=frame.index,
            )
        )

    return frame


def is_valid_bed(
    bed: Union[str, os.PathLike, pd.DataFrame, pr.PyRanges], verbose=True
) -> bool:
    from loguru import logger

    """Return True when the first non-empty row has at least three BED columns."""

    try:
        df = _read_bed_dataframe(bed, nrows=1)
    except FileNotFoundError:
        if verbose:
            logger.warning(f"Bed file: {bed} not found")
        return False
    except pd.errors.EmptyDataError:
        if verbose:
            logger.warning(f"Bed file: {bed} is empty")
        return False
    except Exception as e:
        if verbose:
            logger.warning(f"Exception raised {e}")
        return False

    return df.shape[1] >= 3


def bed_has_name(
    bed: Union[str, os.PathLike, pd.DataFrame, pr.PyRanges]
) -> bool:
    """Return True when the first non-empty row has at least four BED columns."""

    try:
        df = _read_bed_dataframe(bed, nrows=1)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return False

    return df.shape[1] >= 4


def bed_has_duplicate_names(
    bed: Union[str, os.PathLike, pd.DataFrame, pr.PyRanges]
) -> bool:
    """Return True when a BED-like input has duplicate name values."""

    df = convert_bed_to_dataframe(bed)
    if "name" not in df.columns or df.empty:
        return False

    return df["name"].dropna().duplicated().any()


def hash_column(col: Iterable, hash_type=64) -> list:
    """
    Convinience function to perform hashing using xxhash on an iterable.

    Function is **not** vectorised.
    """
    import xxhash

    hash_dict = {
        32: xxhash.xxh32_intdigest,
        64: xxhash.xxh64_intdigest,
        128: xxhash.xxh128_intdigest,
    }

    hash_func = hash_dict.get(hash_type)

    return [hash_func(v) for v in col]


def split_intervals_on_chrom(
    intervals: Union[str, os.PathLike, pd.DataFrame, pr.PyRanges]
) -> dict:
    """Creates dictionary from bed file with the chroms as keys"""

    intervals = convert_bed_to_dataframe(intervals)
    if intervals.empty or "chrom" not in intervals.columns:
        return {}

    return {chrom: df for chrom, df in intervals.groupby("chrom")}


def intersect_bins(
    bins_1: pd.DataFrame, bins_2: pd.DataFrame, **bedtools_kwargs
) -> pd.DataFrame:
    """Intersect two interval tables and return a labeled pandas DataFrame."""

    left = _prepare_intersection_frame(bins_1, name_prefix="region_1")
    right = _prepare_intersection_frame(bins_2, name_prefix="region_2")

    required_columns = [
        "chrom_1",
        "start_1",
        "end_1",
        "name_1",
        "chrom_2",
        "start_2",
        "end_2",
        "name_2",
        "overlap",
    ]

    if left.empty or right.empty:
        return pd.DataFrame(columns=required_columns)

    left = left.rename(
        columns={"chrom": "chrom_1", "start": "start_1", "end": "end_1", "name": "name_1"}
    ).copy()
    right = right.rename(
        columns={"chrom": "chrom_2", "start": "start_2", "end": "end_2", "name": "name_2"}
    ).copy()

    slack = int(bedtools_kwargs.get("slack", 0) or 0)
    if slack:
        left["start_1"] = (left["start_1"] - slack).clip(lower=0)
        left["end_1"] = left["end_1"] + slack
        right["start_2"] = (right["start_2"] - slack).clip(lower=0)
        right["end_2"] = right["end_2"] + slack

    left["_key"] = 1
    right["_key"] = 1
    df_intersect = left.merge(right, on="_key", how="inner").drop(columns="_key")
    df_intersect = df_intersect[df_intersect["chrom_1"] == df_intersect["chrom_2"]]

    if "strandedness" in bedtools_kwargs or bedtools_kwargs.get("s"):
        strandedness = bedtools_kwargs.get("strandedness")
        if bedtools_kwargs.get("s"):
            strandedness = "same"
        if "strand_1" in df_intersect.columns and "strand_2" in df_intersect.columns:
            if strandedness == "same":
                df_intersect = df_intersect[
                    df_intersect["strand_1"] == df_intersect["strand_2"]
                ]
            elif strandedness == "opposite":
                df_intersect = df_intersect[
                    df_intersect["strand_1"] != df_intersect["strand_2"]
                ]

    overlap = (
        df_intersect[["end_1", "end_2"]].min(axis=1)
        - df_intersect[["start_1", "start_2"]].max(axis=1)
    )
    df_intersect = df_intersect.loc[overlap > 0, required_columns[:-1]].copy()
    if df_intersect.empty:
        return pd.DataFrame(columns=required_columns)

    df_intersect["overlap"] = overlap.loc[df_intersect.index]
    return df_intersect[required_columns]


def load_dict(fn, format: str, dtype: str = "int") -> dict:
    """Convinence function to load gziped json/pickle file using xopen."""

    import itertools

    import ujson
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
    import ujson

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
    import time
    from datetime import timedelta

    from loguru import logger

    def wrapper(f):
        @wraps(f)
        def wrapped(*args, **kwargs):
            time_start = time.perf_counter()
            result = f(*args, **kwargs)
            time_end = time.perf_counter()

            time_taken = timedelta(seconds=(time_end - time_start))
            logger.info(f"Completed {task_name} in {time_taken} (hh:mm:ss.ms)")
            return result

        return wrapped

    return wrapper


def convert_to_bedtool(
    bed: Union[str, os.PathLike, pd.DataFrame, pr.PyRanges]
) -> pr.PyRanges:
    """Legacy helper preserved as a PyRanges conversion wrapper."""

    return convert_bed_to_pr(bed)


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
        "norm": "Replicates_Scaled",
        "summary": "Samples_Summarised",
        "subtraction": "Samples_Compared",
    }
    categories = []
    for index, value in ser.iteritems():
        for key in mapping:
            if key in value:
                categories.append(mapping[key])

    return categories


def convert_bed_to_pr(
    bed: Union[
        str,
        os.PathLike,
        pd.DataFrame,
        pr.PyRanges,
    ],
) -> pr.PyRanges:
    """Converts a bed file to a PyRanges object.
    Args:
        bed (Union[str, os.PathLike, pd.DataFrame, pr.PyRanges]): Bed file to convert.
    Returns:
        pr.PyRanges: PyRanges object.
    """

    if isinstance(bed, pr.PyRanges):
        return bed

    df = convert_bed_to_dataframe(bed)
    if df.empty:
        return pr.PyRanges()

    df = _standardize_bed_columns(df, capitalized=True)
    for column in ("Start", "End"):
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="raise")
    if "Name" in df.columns:
        df["Name"] = df["Name"].astype("category")

    return pr.PyRanges(df)


def convert_bed_to_dataframe(
    bed: Union[str, os.PathLike, pd.DataFrame, pr.PyRanges],
    ignore_ray_objrefs=False,
) -> pd.DataFrame:
    """Converts a bed like object (including paths to bed files) to a pd.DataFrame"""
    import ray
    from loguru import logger

    if isinstance(bed, (str, os.PathLike)):
        try:
            bed_conv = _read_bed_dataframe(bed)
        except FileNotFoundError:
            logger.warning(f"File {bed} not found")
            bed_conv = pd.DataFrame()
        except pd.errors.EmptyDataError:
            logger.warning(f"File {bed} is empty")
            bed_conv = pd.DataFrame()

    elif isinstance(bed, pd.DataFrame):
        bed_conv = bed.copy()

    elif isinstance(bed, pr.PyRanges):
        bed_conv = _standardize_bed_columns(bed.copy(), capitalized=False)

    elif isinstance(bed, ray.ObjectRef):
        if ignore_ray_objrefs:
            logger.warning("Assuming ObjectRef is a PyRanges")
            return bed
        else:
            bed = ray.get(bed)
            bed_conv = convert_bed_to_dataframe(bed)

    bed_conv = _standardize_bed_columns(bed_conv, capitalized=False)

    return bed_conv


def is_tabix(file: str):
    import pysam
    from loguru import logger

    _is_tabix = False

    try:
        tbx = pysam.TabixFile(file)
        _chroms = tbx.contigs
        _is_tabix = True

    except OSError as e:
        logger.warn(e)

    return _is_tabix


def format_coordinates(coordinates: Union[str, os.PathLike]) -> pr.PyRanges:
    """Convert coordinates supplied in string format or a BED file to PyRanges.

    Args:
        coordinates (Union[str, os.PathLike]): Coordinates in the form chr:start-end/path.
    Raises:
        ValueError: Inputs must be supplied in the correct format.

    Returns:
        pr.PyRanges: PyRanges object containing the required coordinates.
    """

    coordinates = str(coordinates)
    pattern_genomic_coord = re.compile(
        r"^(chr[0-2xXyYmM][0-9]*):(\d+)-(\d+)(?:\s+(\S+))?$"
    )

    match = pattern_genomic_coord.match(coordinates)
    if match:
        chrom, start, end, name = match.groups()
        if not name:
            name = "region_0"

        return pr.PyRanges(
            pd.DataFrame(
                {
                    "Chromosome": [chrom],
                    "Start": [int(start)],
                    "End": [int(end)],
                    "Name": [name],
                }
            )
        )

    path_name = Path(coordinates).name.lower()
    if path_name.endswith((".bed", ".bed.gz", ".bed.bgz")):
        if is_valid_bed(coordinates):
            bed_df = convert_bed_to_dataframe(coordinates)
            if bed_has_name(bed_df):
                return convert_bed_to_pr(bed_df)

            bed_df = bed_df[["chrom", "start", "end"]].copy()
            bed_df = bed_df.reset_index(drop=True)
            bed_df["name"] = bed_df.index.map(lambda idx: f"region_{idx}")
            return convert_bed_to_pr(bed_df)

        raise ValueError("Invalid bed file supplied.")

    raise ValueError(
        """Provide coordinates in the form chr[NUMBER]:[START]-[END]/BED file"""
    )


def convert_interval_to_coords(
    interval: Union[dict, pd.Series], named=False
) -> Tuple[str, str]:
    """Converts interval object to standard genomic coordinates.

    e.g. chr1:1000-2000

    Args:
        interval (Union[dict, pd.Series]): Interval to convert.

    Returns:
        Tuple: Genomic coordinates in the format chr:start-end
    """
    chrom = interval.get("chrom", interval.get("Chromosome"))
    start = interval.get("start", interval.get("Start"))
    end = interval.get("end", interval.get("End"))
    name = interval.get("name", interval.get("Name", "Unnammed"))

    if not named:
        return ("Unnammed", f"{chrom}:{start}-{end}")
    else:
        return (name, f"{chrom}:{start}-{end}")


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
    from loguru import logger

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
        logger.debug(f"File extension {ext} is not supported")
        raise e


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


def get_restriction_site(restriction_enzyme: str):
    """
    Gets the restriction site for a given restriction enzyme.

    Can be either the name of the restriction enzyme or the restriction site itself.
    The restriction site will just be returned if it is a valid DNA sequence.

    Args:
        restriction_enzyme: Name of restriction enzyme or restriction site.

    Returns:
        Restriction site.

    Raises:
        ValueError: If restriction enzyme is not found.

    """

    if re.match(r"^[ACGTacgt]+$", restriction_enzyme):
        return restriction_enzyme

    import Bio.Restriction

    all_enzymes = {e.lower(): e for e in Bio.Restriction.AllEnzymes.as_string()}
    if restriction_enzyme.lower() not in all_enzymes:
        raise ValueError(f"Restriction enzyme {restriction_enzyme} not found.")
    else:
        return Bio.Restriction.AllEnzymes.get(
            all_enzymes[restriction_enzyme.lower()]
        ).site
