import pathlib
import re
from typing import Dict, List

import pandas as pd
import pyranges as pr
import logging

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


def get_design_matrix(fastqs: List[str, pathlib.Path]):
    
    df = pd.DataFrame(fastqs, columns=["fn"])
    df["filename"] = df["fn"].apply(str).str.split(".fastq").str[0]
    df["samplename"] = df["filename"].str.extract(r".*/(.*?)_R?[12].fastq.*")
    df["condition"] = df["samplename"].str.split(".fastq").str[0].str.split("_").str[-1]

    if df["condition"].isna().any():
        logging.warn("Failed to identify conditions from fastq files. Please format as SAMPLENAME_CONDITION_READ.fastq(.gz)")
        df["condition"].fillna("UNKNOWN")
    
    return df[["sample_name", "condition"]].drop_duplicates()

def get_binsizes(config):
    try:
        binsizes = [
            int(bs)
            for bs in re.split(r"[,;]\s*|\s+", str(config["analysis"].get("bin_sizes", "0")))
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
        import coolbox
    except ImportError as e:
        logging.warning(
            "Plotting capabilities not installed. For plotting please run: pip install capcruncher[plotting]"
        )
        return False
    
    return utils.is_valid_bed(config["plot"].get("coordinates"), verbose=False)
    
