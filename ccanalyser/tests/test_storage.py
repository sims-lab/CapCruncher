import os
import sys
import pandas as pd
import cooler
import pytest
import click
from click.testing import CliRunner

from ccanalyser.tools.storage import GenomicBinner, CoolerBinner, create_cooler_cc


# Pre-run setup
dir_test = os.path.realpath(os.path.dirname(__file__))
dir_package = os.path.dirname(dir_test)
dir_data = os.path.join(dir_package, "data")


def test_make_cooler():
    pixels = pd.read_csv(
        os.path.join(dir_data, "test", "RS411_0uM_1.BCL2.tsv.gz"), sep="\t"
    )
    bins = pd.read_csv(
        os.path.join(dir_data, "test", "genome.digest.bed.gz"),
        sep="\t",
        names=["chrom", "start", "end", "name"],
    )
    oligos = os.path.join(dir_data, "test", "capture_oligos.bed")
    output_prefix = "test/cooler"
    outfile = "test/cooler.BCL2.hdf5"

    if os.path.exists(outfile):
        os.unlink(outfile)

    create_cooler_cc(output_prefix, bins, pixels, "BCL2", oligos)

    assert os.path.exists(outfile)

    clr = cooler.Cooler(outfile)

    assert clr.pixels()[:].shape == pixels.shape
    assert clr.pixels()[:]["count"].sum() == pixels["count"].sum()


def test_binning():

    cooler_fn = "test/cooler.BCL2.hdf5"
    outfile = "test/cooler.BCL2.binned.hdf5"
    cb = CoolerBinner(cooler_fn, binsize=2500, n_cores=8)
    cb.to_cooler(outfile, normalise=False, scale_factor=1e6)

    assert os.path.exists(outfile)
