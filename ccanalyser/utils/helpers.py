import os
import pandas as pd
from typing import Union
import pybedtools
from pybedtools import BedTool
import re


def is_on(param):
    values = ["true", "t", "on", "yes", "y", "1"]
    if str(param).lower() in values:
        return True

def is_off(param):
    values = ['', 'None', 'none', 'F', 'f']
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

    if re.match(r'[GgAaTtCc]+', recognition_site):
        # This matches a DNA sequence so just convert to upper case and return
        return recognition_site.upper()
    elif recognition_site.lower() in known_enzymes:
        return known_enzymes.get(recognition_site.lower())
    else:
        raise ValueError("No restriction site or recognised enzyme provided")

