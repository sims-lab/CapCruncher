import pathlib
import re
from typing import Dict, List, Union

import pandas as pd
import pyranges as pr
import logging
import os

from capcruncher import utils


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


def convert_empty_yaml_entry_to_string(param: str) -> str:
    """
    Converts empty yaml entries to string
    """
    if is_none(param):
        return ""
    else:
        return param


def format_config_dict(config: Dict) -> Dict:
    """
    Formats the config dictionary to ensure that all entries are strings.

    """
    for key, value in config.items():
        config[key] = convert_empty_yaml_entry_to_string(value)

    return config


def get_design_matrix(fastqs: List[Union[str, pathlib.Path]]):

    df = pd.DataFrame(fastqs, columns=["fn"])
    df["filename"] = df["fn"].apply(str).str.split(".fastq").str[0]
    df["sample"] = df["filename"].str.extract(r".*/(.*?)_R?[12].fastq.*")
    df["condition"] = df["sample"].str.split(".fastq").str[0].str.split("_").str[-1]

    if df["condition"].isna().any():
        logging.warn(
            "Failed to identify conditions from fastq files. Please format as sample_CONDITION_READ.fastq(.gz)"
        )
        df["condition"].fillna("UNKNOWN")

    return df[["sample_name", "condition"]].drop_duplicates()


def get_binsizes(config):
    try:
        binsizes = [
            int(bs)
            for bs in re.split(
                r"[,;]\s*|\s+", str(config["analysis"].get("bin_sizes", "0"))
            )
        ]
    except ValueError:
        binsizes = []

    return binsizes


def get_blacklist(config):
    if utils.is_valid_bed(config["analysis_optional"].get("blacklist"), verbose=False):
        blacklist = config["analysis_optional"]["blacklist"]
        return blacklist


def has_high_viewpoint_number(viewpoints: str, config: Dict):
    n_viewpoints = pr.read_bed(viewpoints).df.shape[0]
    if n_viewpoints > 100:
        if not config["analysis_optional"].get("force_bigwig_generation", False):
            return True


def can_perform_plotting(config):
    try:
        pass
    except ImportError:
        logging.warning(
            "Plotting capabilities not installed. For plotting please run: pip install capcruncher[plotting]"
        )
        return False

    return utils.is_valid_bed(config["plot"].get("coordinates"), verbose=False)


def group_files_by_regex(files: List, regex: str):
    df = pd.DataFrame(files, columns=["fn"])
    extracted_substrings = df["fn"].astype(str).str.extract(regex)
    df = df.join(extracted_substrings)
    return (
        df.groupby(extracted_substrings.columns.to_list())
        .agg(list)["fn"]
        .rename("files_grouped")
    )


class GenericFastqSamples:
    def __init__(self, design):

        # Expected columns: sample, fq1, fq2
        self.design = design
        self.design = self.design.assign(
            paired=(~self.design[["fq1", "fq2"]].isna().any(axis=1))
        )

        assert self.design["paired"].all(), "All files must be paired"

    @classmethod
    def from_files(cls, files: List[Union[pathlib.Path, str]]) -> (pd.DataFrame):

        try:
            df = pd.DataFrame(files, columns=["fn"])

            df[["sample", "read"]] = (
                df["fn"].apply(str).str.extract("(?!.*/)?(.*).*_R?([12]).fastq(?!.gz)?")
            )

            df["sample"] = df["sample"].apply(
                lambda p: pathlib.Path(p).name
                if isinstance(p, pathlib.Path)
                else os.path.basename(p)
            )
            df["read"] = "fq" + df["read"]

            df = (
                df.pivot(columns="read", index=["sample"])
                .droplevel(level=0, axis=1)
                .reset_index()
            )

            return cls(design=df)
        except TypeError:
            logging.exception("Samples not named correctly")
            raise ValueError("Could not determine sample name from the file name")

    @property
    def fastq_files(self):
        return sorted([*self.design["fq1"], *self.design["fq2"]])

    @property
    def sample_names_all(self):
        return self.design["sample"].to_list()

    @property
    def experimental_design(self):
        return self.design.assign(
            condition=self.design["sample"].split("_").str[1]
        ).loc[:, ["sample", "condition"]]
