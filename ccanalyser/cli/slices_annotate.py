import itertools

import click
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from pybedtools import BedTool

from ccanalyser.cli import cli
from ccanalyser.utils import bed_has_name, is_valid_bed


def find_intersections(
    a,
    b,
    frac=1e-9,
    column_name="count",
    method="get",
    method_na_values={"get": ".", "count": 0},
):

    if isinstance(a, str):
        a = BedTool(a)
    elif isinstance(a, pd.DataFrame):
        a = BedTool.from_dataframe(a)

    if is_valid_bed(b):

        if method == "get":
            bt_intersections = a.intersect(b, loj=True, f=frac)
        elif method == "count":
            bt_intersections = a.intersect(b, loj=True, c=True, f=frac)
        else:
            raise ValueError("method argument must be in [get|count]")

        return format_intersections(
            intersections=bt_intersections, column_name=column_name
        )

    else:
        return format_intersections(
            intersections=a,
            column_name=column_name,
            failed=True,
            na_value=method_na_values.get(method, np.nan),
        )


def format_intersections(intersections, column_name, failed=False, na_value=0):

    df_intersections = intersections.to_dataframe()

    if not failed:
        # Rename last column to column name
        df_intersections = df_intersections.rename(
            columns={df_intersections.columns[-1]: column_name}
        )

        # Extract only needed columns
        df_intersections = df_intersections[["name", "chrom", "start", column_name]]

    else:
        df_intersections[column_name] = na_value
        df_intersections = df_intersections[["name", "chrom", "start", column_name]]

    return df_intersections.set_index("name")[column_name]


def cycle_argument(arg):
    """Allows for the same argument to be stated once but repeated for all files"""

    if len(arg) == 1:
        arg_cp = arg.copy()
        return itertools.cycle(arg_cp)
    else:
        return arg


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
def slices_annotate(
    actions=None,
    slices=None,
    bed_files=None,
    names=None,
    overlap_fractions=None,
    output=None,
    duplicates="remove",
):

    """Annotates a bam file (converted to bed format) with other bed files"""


    #TODO: Increase speed by splitting files on chromosome and intersecting in parallel

    # Verify bed integrity
    assert is_valid_bed(slices), f"bed - {slices} is invalid"
    assert bed_has_name(slices), f"bed - {slices} does not have a name column"

    print("Formating  bed file")
    # Make base dataframe
    df_bed = df_bed = BedTool(slices).to_dataframe()

    print("Dealing with any duplicate names")
    # Deal with duplicates -- currently very harsh
    if duplicates == "remove":
        df_bed = df_bed.sort_values(["score"], ascending=False).drop_duplicates(
            subset="name", keep="first"
        )

        slices = BedTool.from_dataframe(df_bed[["chrom", "start", "end", "name"]])

    # Run the intersection
    n_actions = len(actions)
    dframes = Parallel(n_jobs=n_actions)(
        delayed(find_intersections)(
            a=slices, b=b2, frac=f, column_name=name, method=action
        )
        for b2, f, name, action in zip(
            bed_files,
            cycle_argument(overlap_fractions),
            names,
            actions,
        )
    )

    # Merge dataframe with annotations
    df_annotations = (
        df_bed.set_index(["name"])
        .join(dframes, how="left")
        .reset_index()
        .rename(columns={"name": "slice_name"})
        .drop(columns=["score", "strand"], errors="ignore")
    )

    # Export to csv
    df_annotations.to_csv(output, sep="\t", index=False)
