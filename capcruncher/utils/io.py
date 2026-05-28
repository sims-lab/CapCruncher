"""File I/O utilities: dict serialisation, cooler URIs, tabix, file-type detection."""

from __future__ import annotations

import os
import pickle
from collections.abc import Iterable
from typing import Any

import pandas as pd

from capcruncher.types import (
    DictDType,
    DictFormat,
    VALID_DICT_DTYPES,
    VALID_DICT_FORMATS,
    validate_choice,
)


def read_dataframes(filenames: Iterable, **kwargs) -> list[pd.DataFrame]:
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

    import json

    from xopen import xopen

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


def get_file_type(fn: os.PathLike) -> str:
    """Determines file type based on extension."""
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
) -> str:
    import re

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


def is_tabix(file: str) -> bool:
    import pysam
    from loguru import logger

    _is_tabix = False

    try:
        tbx = pysam.TabixFile(file)
        _chroms = tbx.contigs
        _is_tabix = True

    except OSError as e:
        logger.warning(e)

    return _is_tabix
