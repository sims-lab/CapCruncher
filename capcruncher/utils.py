import itertools
import os
import pickle
import re
from functools import wraps
from typing import Callable, Iterable, Mapping, Sequence, Tuple, Union

import pandas as pd
import pyranges1 as pr


BED_COLUMN_NAMES = ("chrom", "start", "end", "name", "score", "strand")
PYRANGES_TO_BED_COLUMNS = {
    "Chromosome": "chrom",
    "Start": "start",
    "End": "end",
    "Name": "name",
    "Score": "score",
    "Strand": "strand",
}
BED_TO_PYRANGES_COLUMNS = {value: key for key, value in PYRANGES_TO_BED_COLUMNS.items()}


def _pyranges_to_dataframe(gr: pr.PyRanges) -> pd.DataFrame:
    return pd.DataFrame(gr).copy()


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


def _bed_column_names(n_columns: int) -> list[str]:
    base_columns = list(BED_COLUMN_NAMES[: min(n_columns, len(BED_COLUMN_NAMES))])
    extra_columns = [f"column_{index}" for index in range(len(base_columns) + 1, n_columns + 1)]
    return [*base_columns, *extra_columns]


def _coerce_bedframe(df: pd.DataFrame) -> pd.DataFrame:
    bed = df.copy()

    if all(isinstance(column, int) for column in bed.columns):
        bed.columns = _bed_column_names(bed.shape[1])
    else:
        bed = bed.rename(columns=PYRANGES_TO_BED_COLUMNS)

    required_columns = [column for column in BED_COLUMN_NAMES[:3] if column in bed.columns]
    if len(required_columns) < 3:
        raise IndexError("Wrong number of fields detected check separator or number of columns")

    bed["start"] = pd.to_numeric(bed["start"], errors="raise").astype(int)
    bed["end"] = pd.to_numeric(bed["end"], errors="raise").astype(int)

    return bed


def _read_bed_file(path: Union[str, os.PathLike]) -> pd.DataFrame:
    bed = pd.read_csv(path, sep="\t", header=None, comment="#")
    if bed.empty:
        raise pd.errors.EmptyDataError(f"{path} is empty")

    bed.columns = _bed_column_names(bed.shape[1])
    return _coerce_bedframe(bed)


def _ensure_name_column(df: pd.DataFrame, prefix: str = "region") -> pd.DataFrame:
    bed = df.copy()
    if "name" not in bed.columns:
        bed["name"] = [f"{prefix}_{index}" for index in range(len(bed))]
    else:
        missing_names = bed["name"].isna() | (bed["name"].astype(str).str.len() == 0)
        if missing_names.any():
            bed.loc[missing_names, "name"] = [
                f"{prefix}_{index}" for index in bed.index[missing_names]
            ]

    return bed


def _parse_region(region: Union[str, Sequence[Union[str, int]]]) -> tuple[str, int, int]:
    if isinstance(region, str):
        chrom, coordinates = region.split(":", maxsplit=1)
        start, end = coordinates.split("-", maxsplit=1)
        return chrom, int(start), int(end)

    chrom, start, end = region[:3]
    return str(chrom), int(start), int(end)


def is_valid_bed(bed: Union[str, os.PathLike, pd.DataFrame, pr.PyRanges], verbose=True) -> bool:
    from loguru import logger

    """Returns true if bed file can be opened and has at least 3 columns"""
    try:
        bed = convert_bed_to_dataframe(bed)
        if {"chrom", "start", "end"}.issubset(bed.columns):
            return True

    except Exception as e:
        if not verbose:
            return False

        if isinstance(e, FileNotFoundError):
            logger.warning(f"Bed file: {bed} not found")

        elif isinstance(e, (IndexError, ValueError, pd.errors.ParserError)):
            logger.warning(
                "Wrong number of fields detected check separator or number of columns"
            )

        else:
            logger.warning(f"Exception raised {e}")

    return False


def bed_has_name(bed: Union[str, os.PathLike, pd.DataFrame, pr.PyRanges]) -> bool:
    """Returns true if bed file has at least 4 columns"""
    bed_df = convert_bed_to_dataframe(bed)
    return "name" in bed_df.columns


def bed_has_duplicate_names(bed: Union[str, os.PathLike, pd.DataFrame, pr.PyRanges]) -> bool:
    """Returns true if the bed file contains duplicated names."""
    bed_df = convert_bed_to_dataframe(bed)
    return "name" in bed_df.columns and bed_df["name"].duplicated().any()


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
    return {chrom: df for chrom, df in intervals.groupby("chrom")}


def intersect_bins(
    bins_1: pd.DataFrame, bins_2: pd.DataFrame, **bedtools_kwargs
) -> pd.DataFrame:
    """Intersects two sets of genomic intervals using PyRanges joins."""

    _ = bedtools_kwargs

    pr_1 = convert_bed_to_pr(_ensure_name_column(bins_1))
    pr_2 = convert_bed_to_pr(_ensure_name_column(bins_2))
    df_intersect = _pyranges_to_dataframe(
        pr_1.join_overlaps(pr_2, report_overlap_column="Overlap")
    )

    if df_intersect.empty:
        return pd.DataFrame(
            columns=[
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
        )

    return df_intersect.rename(
        columns={
            "Chromosome": "chrom_1",
            "Start": "start_1",
            "End": "end_1",
            "Name": "name_1",
            "Chromosome_b": "chrom_2",
            "Start_b": "start_2",
            "End_b": "end_2",
            "Name_b": "name_2",
            "Overlap": "overlap",
        }
    )[
        [
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
    ]


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
) -> pd.DataFrame:
    """Returns a dataframe representation of a bed-like object."""
    return convert_bed_to_dataframe(bed)


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
        bed (Union[str, pybedtools.BedTool, pd.DataFrame, pr.PyRanges]): Bed file to convert.
    Returns:
        pr.PyRanges: PyRanges object.
    """

    if isinstance(bed, pr.PyRanges):
        converted = bed

    else:
        try:
            bed_df = convert_bed_to_dataframe(bed)
        except (FileNotFoundError, pd.errors.EmptyDataError):
            from loguru import logger

            logger.warning(f"File {bed} not found")
            return pr.PyRanges()

        converted = (
            _ensure_name_column(bed_df)
            .rename(columns=BED_TO_PYRANGES_COLUMNS)
            .assign(Name=lambda df: df["Name"].astype("category"))
            .pipe(pr.PyRanges)
        )

    return converted


def convert_bed_to_dataframe(
    bed: Union[
        str, os.PathLike, pd.DataFrame, pr.PyRanges
    ],  # noqa: F821
    ignore_ray_objrefs=False,
) -> pd.DataFrame:
    """Converts a bed like object (including paths to bed files) to a pd.DataFrame"""
    import ray
    from loguru import logger

    if isinstance(bed, (str, os.PathLike)):
        bed_conv = _read_bed_file(bed)

    elif isinstance(bed, pd.DataFrame):
        bed_conv = _coerce_bedframe(bed)

    elif isinstance(bed, pr.PyRanges):
        bed_conv = _coerce_bedframe(_pyranges_to_dataframe(bed))

    elif isinstance(bed, ray.ObjectRef):
        if ignore_ray_objrefs:
            logger.warning("Assuming ObjectRef is a PyRanges")
            bed_conv = bed
        else:
            bed = ray.get(bed)
            bed_conv = convert_bed_to_dataframe(bed)

    return bed_conv


def bam_to_bed_dataframe(bam: Union[str, os.PathLike]) -> pd.DataFrame:
    import pysam

    rows = []
    with pysam.AlignmentFile(bam, "rb") as alignment_file:
        for alignment in alignment_file.fetch(until_eof=True):
            if alignment.is_unmapped:
                continue

            rows.append(
                {
                    "chrom": alignment.reference_name,
                    "start": alignment.reference_start,
                    "end": alignment.reference_end,
                    "name": alignment.query_name,
                    "score": alignment.mapping_quality,
                    "strand": "-" if alignment.is_reverse else "+",
                }
            )

    return pd.DataFrame(rows, columns=list(BED_COLUMN_NAMES))


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


def fetch_bed_intervals(
    bed: Union[str, os.PathLike],
    region: Union[str, Sequence[Union[str, int]]],
) -> pd.DataFrame:
    import pysam

    chrom, start, end = _parse_region(region)

    if is_tabix(bed):
        with pysam.TabixFile(bed) as tabix_file:
            rows = list(tabix_file.fetch(chrom, start, end, parser=pysam.asTuple()))

        if not rows:
            return pd.DataFrame(columns=list(BED_COLUMN_NAMES[:4]))

        fetched = pd.DataFrame(rows, columns=_bed_column_names(len(rows[0])))
        fetched = _coerce_bedframe(fetched)
    else:
        fetched = convert_bed_to_dataframe(bed)
        fetched = fetched.loc[
            lambda df: (df["chrom"] == chrom) & (df["start"] < end) & (df["end"] > start)
        ]

    return fetched.reset_index(drop=True)


def format_coordinates(coordinates: Union[str, os.PathLike]) -> pd.DataFrame:
    """Converts coordinates supplied in string format or a .bed file to a dataframe.

    Args:
        coordinates (Union[str, os.PathLike]): Coordinates in the form chr:start-end/path.
    Raises:
        ValueError: Inputs must be supplied in the correct format.

    Returns:
        pd.DataFrame: DataFrame containing the required coordinates.
    """

    coordinates = str(coordinates)
    pattern_genomic_coord = re.compile(r"chr[0-2xXyYmM][0-9]*:\d+-\d+(\s\w)*$")
    pattern_bed_file = re.compile(r"(.*).bed")

    if pattern_genomic_coord.match(coordinates):
        coordinates_split = re.split(":|-", coordinates)
        if len(coordinates_split) < 4:
            coordinates_split.append("region_0")

        bed = pd.DataFrame(
            [coordinates_split], columns=["chrom", "start", "end", "name"]
        ).assign(start=lambda df: df["start"].astype(int), end=lambda df: df["end"].astype(int))

    elif pattern_bed_file.match(coordinates):
        if is_valid_bed(coordinates):
            bed = convert_bed_to_dataframe(coordinates)
            bed = _ensure_name_column(bed.reset_index(drop=True))[
                ["chrom", "start", "end", "name"]
            ]
        else:
            raise ValueError("Invalid bed file supplied.")

    else:
        raise ValueError(
            """Provide coordinates in the form chr[NUMBER]:[START]-[END]/BED file"""
        )

    return bed


def convert_interval_to_coords(
    interval: Union[Mapping, pd.Series], named=False
) -> Tuple[str]:
    """Converts interval object to standard genomic coordinates.

    e.g. chr1:1000-2000

    Args:
        interval (Union[Mapping, pd.Series]): Interval to convert.

    Returns:
        Tuple: Genomic coordinates in the format chr:start-end
    """
    if hasattr(interval, "_asdict"):
        interval = interval._asdict()

    chrom = interval.get("chrom", interval.get("Chromosome"))
    start = interval.get("start", interval.get("Start"))
    end = interval.get("end", interval.get("End"))
    name = interval.get("name", interval.get("Name", "Unnammed"))

    if not named:
        return (
            "Unnammed",
            f"{chrom}:{start}-{end}",
        )
    else:
        return (
            name,
            f"{chrom}:{start}-{end}",
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
