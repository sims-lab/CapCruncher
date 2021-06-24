import itertools
import sys
import warnings
from typing import Tuple, Union

warnings.simplefilter("ignore")
import os

import pandas as pd
from capcruncher.tools.annotate import BedIntersection
from capcruncher.utils import (bed_has_name, convert_bed_to_dataframe,
                              is_valid_bed)
from pybedtools import BedTool


def cycle_argument(arg):
    """Allows for the same argument to be stated once but repeated for all files"""

    if len(arg) == 1:
        return itertools.cycle((arg[0],))
    else:
        return arg


def remove_duplicates_from_bed(bed: Union[str, BedTool, pd.DataFrame]) -> BedTool:
    """
    Simple removal of duplicated entries from bed file.

    If a "score" field is present a higher scored entry is prioritised.

    Args:
     bed (Union[str, BedTool, pd.DataFrame]): Bed object to deduplicate

    Returns:
     BedTool: BedTool with deduplicated names
    """

    df = convert_bed_to_dataframe(bed).sample(frac=1)

    if "score" in df.columns:
        df = df.sort_values(["score"], ascending=False)

    return (
        df.drop_duplicates(subset="name", keep="first")
        .sort_values(["chrom", "start"])[["chrom", "start", "end", "name"]]
        .pipe(BedTool.from_dataframe)
    )



def annotate(
    slices: os.PathLike,
    actions: Tuple = None,
    bed_files: Tuple = None,
    names: Tuple = None,
    overlap_fractions: Tuple = None,
    output: os.PathLike = None,
    duplicates: str = "remove",
    n_cores: int = 8,
    invalid_bed_action: str = "error",
):
    """
    Annotates a bed file with other bed files using bedtools intersect.

    Whilst bedtools intersect allows for interval names and counts to be used for annotating intervals, this command
    provides the ability to annotate intervals with both interval names and counts at the same time. As the pipeline allows
    for empty bed files, this command has built in support to deal with blank/malformed bed files and will return default N/A values.

    Prior to interval annotation, the bed file to be intersected is validated and duplicate entries/multimapping reads are removed 
    to ensure consistent annotations and prevent issues with reporter identification. 

    \f
    Args:
     slices (os.PathLike): Input bed file.
     actions (Tuple, optional): Methods to use for annotation. Choose from (get|count). Defaults to None.
     bed_files (Tuple, optional): Bed files to intersect with the bed file to be annotated. Defaults to None.
     names (Tuple, optional): Column names for output tsv file. Defaults to None.
     overlap_fractions (Tuple, optional): Minimum overlap fractions required to call an intersection. Defaults to None.
     output (os.PathLike, optional): Output file path for annotated .tsv file. Defaults to None.
     duplicates (str, optional): Method to deal with multimapping reads/duplicate bed names.
                                 Currently, "remove" is the only supported option. Defaults to "remove".
     n_cores (int, optional): Number of corese to use for intersection. Bed files are split by chromosome for faster intersection.
                              Defaults to 8.
     invalid_bed_action (str, optional): Action to deal with invalid bed files. Choose from (ignore|error) .These can be ignored by setting to "ignore". Defaults to 'error'.

    Raises:
     NotImplementedError: Only supported option for duplicate bed names is remove.
    """

    # If reading from stdin
    if slices == '-':
        slices = pd.read_csv(sys.stdin, sep='\t', header=None).pipe(BedTool.from_dataframe)

    print('Validating bed file')

    # Check if valid bed format
    if not is_valid_bed(slices):
        raise ValueError(f"bed - {slices} is invalid")
    
    # Check if name column present
    if not bed_has_name(slices):
        raise ValueError(f"bed - {slices} does not have a name column")
    
    # Check if the right number of items present
    if not len(names) == len(bed_files) == len(actions):
        raise IndexError("Wrong number of column names/files/actions provided, check command")

    print('Dealing with duplicates in the bed file')

    # Deal with multimapping reads.
    if duplicates == "remove":  
        slices = remove_duplicates_from_bed(slices)
    else:
        raise NotImplementedError(
            "Only supported option at present is to remove duplicates"
        )

    print('Performing intersection')

    # Perform intersections
    intersection_series = []
    for bed, name, action, fraction in zip(
        bed_files, names, actions, cycle_argument(overlap_fractions)
    ):

        bi = BedIntersection(
            bed1=slices,
            bed2=bed,
            intersection_name=name,
            intersection_method=action,
            intersection_min_frac=fraction,
            intersection_split_chrom=True if n_cores > 1 else False,
            n_cores=n_cores,
            invalid_bed_action=invalid_bed_action,
        )

        intersection_series.append(bi.intersection)

    print('Merging annotations')
    # Merge intersections with slices
    df_annotation = (
        convert_bed_to_dataframe(slices)
        .set_index("name")
        .join(intersection_series)
        .reset_index()
        .rename(columns={"name": "slice_name"})
    )

    print('Writing annotations to file.')
    # Export to tsv
    df_annotation.to_csv(output, sep="\t", index=False)
