import itertools
import os
import pickle
import re
from pathlib import Path
from functools import wraps
from typing import Callable, Iterable

import pandas as pd
import pyranges1 as pr
from capcruncher.types import (
    DictDType,
    DictFormat,
    VALID_DICT_DTYPES,
    VALID_DICT_FORMATS,
    validate_choice,
)

type BedInput = str | os.PathLike | pd.DataFrame | pr.PyRanges

BED_COLUMN_NAMES = [
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

BED_COLUMN_CASE = {
    "chrom": "Chromosome",
    "start": "Start",
    "end": "End",
    "name": "Name",
    "score": "Score",
    "strand": "Strand",
    "thick_start": "ThickStart",
    "thick_end": "ThickEnd",
    "item_rgb": "ItemRGB",
    "block_count": "BlockCount",
    "block_sizes": "BlockSizes",
    "block_starts": "BlockStarts",
}

BED_COLUMN_ALIASES = {
    "chrom": "chrom",
    "chromosome": "chrom",
    "start": "start",
    "end": "end",
    "name": "name",
    "score": "score",
    "strand": "strand",
    "thickstart": "thick_start",
    "thickend": "thick_end",
    "itemrgb": "item_rgb",
    "blockcount": "block_count",
    "blocksizes": "block_sizes",
    "blockstarts": "block_starts",
}

INTERSECT_COLUMNS = [
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
            continue

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
        return f"{bp}bp"
    if (bp / 1e3) < 1000:
        return f"{bp / 1e3}kb"
    return f"{bp / 1e6}mb"


def _read_bed_dataframe(bed: BedInput, nrows=None) -> pd.DataFrame:
    if isinstance(bed, pr.PyRanges):
        return bed.copy()

    if isinstance(bed, pd.DataFrame):
        return bed.copy()

    df = pd.read_csv(
        bed,
        sep="\t",
        header=None,
        comment="#",
        nrows=nrows,
    )

    df.columns = [
        BED_COLUMN_NAMES[i] if i < len(BED_COLUMN_NAMES) else f"col_{i}"
        for i in range(df.shape[1])
    ]
    return df


def _standardize_bed_columns(
    df: pd.DataFrame, capitalized: bool = False
) -> pd.DataFrame:
    rename_map = {}
    for column in df.columns:
        alias_key = re.sub(r"[^a-z0-9]", "", str(column).lower())
        canonical = BED_COLUMN_ALIASES.get(alias_key)
        if canonical:
            rename_map[column] = BED_COLUMN_CASE[canonical] if capitalized else canonical

    return df.rename(columns=rename_map)


def _prepare_intersection_frame(
    df: BedInput, name_prefix: str
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


def is_valid_bed(bed: BedInput, verbose=True) -> bool:
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


def bed_has_name(bed: BedInput) -> bool:
    """Return True when the first non-empty row has at least four BED columns."""

    try:
        df = _read_bed_dataframe(bed, nrows=1)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return False

    return df.shape[1] >= 4


def bed_has_duplicate_names(bed: BedInput) -> bool:
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
    if hash_func is None:
        raise ValueError(f"Unsupported hash type: {hash_type}")

    return [hash_func(v) for v in col]


def split_intervals_on_chrom(intervals: BedInput) -> dict:
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

    if left.empty or right.empty:
        return pd.DataFrame(columns=INTERSECT_COLUMNS)

    left = left.copy()
    right = right.copy()

    slack = int(bedtools_kwargs.get("slack", 0) or 0)
    if slack:
        left["start"] = (left["start"] - slack).clip(lower=0)
        left["end"] = left["end"] + slack
        right["start"] = (right["start"] - slack).clip(lower=0)
        right["end"] = right["end"] + slack

    strandedness = bedtools_kwargs.get("strandedness")
    if bedtools_kwargs.get("s"):
        strandedness = "same"

    joined = convert_bed_to_pr(left).join_overlaps(
        convert_bed_to_pr(right),
        strand_behavior="ignore",
        suffix="_2",
        report_overlap_column="overlap",
    )
    df_intersect = joined.copy()
    if df_intersect.empty:
        return pd.DataFrame(columns=INTERSECT_COLUMNS)

    if strandedness in {"same", "opposite"} and {"Strand", "Strand_2"}.issubset(
        df_intersect.columns
    ):
        if strandedness == "same":
            df_intersect = df_intersect[
                df_intersect["Strand"] == df_intersect["Strand_2"]
            ]
        else:
            df_intersect = df_intersect[
                df_intersect["Strand"] != df_intersect["Strand_2"]
            ]
        if df_intersect.empty:
            return pd.DataFrame(columns=INTERSECT_COLUMNS)

    return pd.DataFrame(
        {
            "chrom_1": df_intersect["Chromosome"],
            "start_1": df_intersect["Start"],
            "end_1": df_intersect["End"],
            "name_1": df_intersect["Name"],
            "chrom_2": df_intersect["Chromosome"],
            "start_2": df_intersect["Start_2"],
            "end_2": df_intersect["End_2"],
            "name_2": df_intersect["Name_2"],
            "overlap": df_intersect["overlap"],
        },
        columns=INTERSECT_COLUMNS,
    )


def load_dict(
    fn: os.PathLike,
    format: DictFormat | str = DictFormat.JSON,
    dtype: DictDType | str = DictDType.INT,
) -> dict | set:
    """Load a gzipped JSON or pickle mapping with validated key/value dtype conversion."""

    import itertools

    import json
    from xopen import xopen

    format = validate_choice(format, VALID_DICT_FORMATS, "format")
    dtype = validate_choice(dtype, VALID_DICT_DTYPES, "dtype")

    d: dict | set
    if format == DictFormat.JSON:
        with xopen(fn) as r:
            d = json.load(r)
    elif format == DictFormat.PICKLE:
        with xopen(fn, "rb") as r:
            d = pickle.load(r)
    else:
        raise ValueError(f"Unsupported dictionary format: {format}")

    key_sample = list(itertools.islice(d, 50))
    dtype_converters = {
        DictDType.INT: int,
        DictDType.STR: str,
    }
    required_dtype = dtype_converters[dtype]

    if all(isinstance(k, required_dtype) for k in key_sample):
        return d
    if isinstance(d, set):
        return {required_dtype(k) for k in d}
    if isinstance(d, dict):
        return {
            required_dtype(k): required_dtype(v) if v else None for k, v in d.items()
        }
    raise TypeError(f"Unsupported serialized object type: {type(d)!r}")


def save_dict(
    obj: dict | set, fn: os.PathLike, format: DictFormat | str = DictFormat.JSON
) -> os.PathLike:
    """Save a dictionary or set as gzipped JSON or pickle."""

    from xopen import xopen
    import json

    format = validate_choice(format, VALID_DICT_FORMATS, "format")

    if format == DictFormat.JSON:
        with xopen(fn, "w") as w:
            if isinstance(obj, set):
                d = dict.fromkeys(obj)
            else:
                d = obj
            json.dump(d, w)
    elif format == DictFormat.PICKLE:
        with xopen(fn, "wb") as w:
            pickle.dump(obj, w)
    else:
        raise ValueError(f"Unsupported dictionary format: {format}")

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


def convert_bed_to_pr(bed: BedInput) -> pr.PyRanges:
    """Convert a BED-like object to a PyRanges object.

    Args:
        bed: BED path, pandas DataFrame, or PyRanges object.

    Returns:
        PyRanges object.
    """

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


def convert_bed_to_dataframe(bed: BedInput) -> pd.DataFrame:
    """Converts a BED-like object to a DataFrame-style interval table.

    PyRanges1 frames are pandas DataFrame subclasses, so in-memory PyRanges
    inputs are copied directly and manipulated with pandas methods.
    """
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

    elif isinstance(bed, pr.PyRanges):
        bed_conv = bed.copy()

    elif isinstance(bed, pd.DataFrame):
        bed_conv = bed.copy()

    else:
        raise TypeError(f"Unsupported BED input type: {type(bed)!r}")

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


def format_coordinates(coordinates: str | os.PathLike) -> pr.PyRanges:
    """Convert coordinates supplied in string format or a BED file to PyRanges.

    Args:
        coordinates: Coordinates in the form chr:start-end or a BED path.
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
    interval: dict | pd.Series, named: bool = False
) -> tuple[str, str]:
    """Converts interval object to standard genomic coordinates.

    e.g. chr1:1000-2000

    Args:
        interval: Interval to convert.

    Returns:
        Pair of name and genomic coordinates in the format chr:start-end.
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


def get_cooler_uri(
    store: os.PathLike | str, viewpoint: str, resolution: str | int | None
):
    store = os.fspath(store)
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
