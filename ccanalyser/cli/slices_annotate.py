import itertools
from typing import Union
import warnings
warnings.simplefilter('ignore')

import click
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from pybedtools import BedTool

from ccanalyser.cli import cli
from ccanalyser.tools.annotate import BedIntersection
from ccanalyser.utils import bed_has_name, convert_bed_to_dataframe, convert_to_bedtool, is_valid_bed
from natsort import natsorted

def cycle_argument(arg):
    """Allows for the same argument to be stated once but repeated for all files"""

    if len(arg) == 1:
        return itertools.cycle((arg[0],))
    else:
        return arg


def remove_duplicates_from_bed(bed):

    if isinstance(bed, str):
        df = BedTool(bed).to_dataframe()
    elif isinstance(bed, BedTool):
        df = bed.to_dataframe()
    
    if 'score' in df.columns:
        df = df.sort_values(["score"], ascending=False)


    return (df.drop_duplicates(subset="name", keep="first")
              .sort_values(['chrom', 'start'])
              [['chrom', 'start', 'end', 'name']]
              .pipe(BedTool.from_dataframe))


@cli.command()
@click.argument("slices")
@click.option(
    "-a",
    "--actions",
    help="Actions to perform for each bed_files file",
    multiple=True,
    type=click.Choice(
        ["get", "count"],
    ),
)
@click.option(
    "-b", "--bed_files", help="Bed files to intersect with slices", multiple=True
)
@click.option("-n", "--names", help="Names to use as column names", multiple=True)
@click.option(
    "-f",
    "--overlap_fractions",
    help="Minimum overlap fractions",
    multiple=True,
    default=[
        1e-9,
    ],
    type=click.FLOAT,
)
@click.option(
    "-o",
    "--output",
    help="Output file name",
    default="annotated.slices.tsv.gz",
)
@click.option(
    "--duplicates",
    help="Method to use for reconciling duplicate (i.e. multimapping) slices",
    type=click.Choice(["remove"]),
    default="remove",
)
@click.option(
    "-p",
    '--n_cores',
    help="Number of cores to use for intersections",
    default=8
)
@click.option(
    "--invalid_bed_action",
    help="Method to deal with invalid bed files",
    default='error',
    type=click.Choice(["ignore", 'error']),
)
def slices_annotate(
    slices,
    actions=None,
    bed_files=None,
    names=None,
    overlap_fractions=None,
    output=None,
    duplicates="remove",
    n_cores=8,
    invalid_bed_action='error'
):

    """Annotates a bam file (converted to bed format) with other bed files"""

    assert is_valid_bed(slices), f"bed - {slices} is invalid" # Make sure file exist and has the correct number of fields
    assert bed_has_name(slices), f"bed - {slices} does not have a name column" # Make sure name column present
    assert len(names) == len(bed_files) == len(actions), 'Wrong number of column names/files/actions provided, check command' # All args present


    if duplicates == "remove": # Deal with multimapping reads (default)
        slices = remove_duplicates_from_bed(slices)
    else:
        raise NotImplementedError('Only supported option at present is to remove duplicates')
   

    # Perform intersections
    intersection_series = []
    for bed, name, action, fraction in zip(bed_files, names, actions, cycle_argument(overlap_fractions)):

        bi = BedIntersection(bed1=slices,
                             bed2=bed,
                             intersection_name=name,
                             intersection_method=action,
                             intersection_min_frac=fraction,
                             n_cores=n_cores,
                             invalid_bed_action=invalid_bed_action)
        
        intersection_series.append(bi.intersection)
    
    

    # Merge intersections with slices
    df_annotation = (convert_bed_to_dataframe(slices)
                           .set_index('name')
                           .join(intersection_series)
                           .reset_index()
                           .rename(columns={'name': 'slice_name'}))


    # Export to csv
    df_annotation.to_csv(output, sep="\t", index=False)
