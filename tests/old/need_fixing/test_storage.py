# import os
# import sys
# import pandas as pd
# import cooler
# import pytest
# import click
# from click.testing import CliRunner
# import glob

# from capcruncher.tools.storage import GenomicBinner, CoolerBinner, create_cooler_cc
# from capcruncher.cli import cli, reporters_store


# # Pre-run setup

# dir_test = os.path.realpath(os.path.dirname(__file__))
# dir_package = os.path.dirname(dir_test)
# dir_data = os.path.join(dir_package, "data")


# @pytest.fixture(scope="module", autouse=True)
# def cleanup():
    
#     yield

#     # Remove previous test files
#     for fn in glob.glob(os.path.join(dir_test, "test", "*.hdf5")):
#         os.unlink(fn)
    
#     if os.path.exists(os.path.join(dir_test, "test", 'genome_binner.pkl')):
#         os.unlink(os.path.join(dir_test, "test", 'genome_binner.pkl'))



# def test_make_cooler():
#     pixels = pd.read_csv(
#       os.path.join(dir_data, "test", "Slc25A37_reporter_counts.tsv.gz"), sep="\t"

#     )
#     bins = pd.read_csv(
#         os.path.join(dir_data, "test", "genome.digest.bed.gz"),
#         sep="\t",
#         names=["chrom", "start", "end", "name"],
#     )
#     oligos = os.path.join(dir_data, "test", "mm9_capture_oligos.bed")
#     output_prefix = os.path.join(dir_test, "test/cooler")
#     outfile = os.path.join(dir_test, "test/cooler.Slc25A37.fragments.hdf5")

#     if os.path.exists(outfile):
#         os.unlink(outfile)

#     create_cooler_cc(output_prefix, bins, pixels, "Slc25A37", oligos, suffix='fragments')

#     assert os.path.exists(outfile)

#     clr = cooler.Cooler(outfile)

#     assert clr.pixels()[:].shape == pixels.shape
#     assert clr.pixels()[:]["count"].sum() == pixels["count"].sum()


# def test_binning():
  
#     cooler_fn = os.path.join(dir_test, "test", "cooler.Slc25A37.fragments.hdf5")
#     outfile = os.path.join(dir_test, "test", "cooler.binned")
#     cb = CoolerBinner(cooler_fn, binsize=2500, n_cores=8)
#     cb.normalise(scale_factor=1e6)
#     cb.to_cooler(outfile)
#     assert os.path.exists(f'{outfile}.Slc25A37.2500.hdf5')


# def test_make_bin_conversion():

#     from capcruncher.tools.storage import GenomicBinner
#     import pickle

#     frags = pd.read_csv(
#         os.path.join(dir_data, "test", "genome.digest.bed.gz"),
#         sep="\t",
#         names=["chrom", "start", "end", "name"],
#     )

#     gb = GenomicBinner(
#             chromsizes=os.path.join(dir_data, "test", "genome.fa.fai"), fragments=frags, binsize=100000)
    
#     gb.bin_conversion_table

#     binners_dict = dict()
#     binners_dict[1000000] = gb


#     with open(os.path.join(dir_test, "test", "genome_binner.pkl"), "wb") as w:
#         pickle.dump(binners_dict, w)
    
# def test_merging():

#     clr_1 = os.path.join(dir_test, "test", "cooler.binned.Slc25A37.2500.hdf5")
#     clr_2 = os.path.join(dir_test, "test/cooler.Slc25A37.fragments.hdf5")

#     import capcruncher.cli.reporters_store

#     outfile = os.path.join(dir_test, "test", "merged.hdf5")
#     reporters_store.merge([clr_1, clr_2], outfile)

#     assert os.path.exists(outfile)



# def test_binning_with_conversion_table():

#     import pickle

#     cooler_fn = os.path.join(dir_test, "test", "cooler.Slc25A37.fragments.hdf5")
#     outfile = os.path.join(dir_test, "test", "cooler.gb.binned")

#     with open(os.path.join(dir_test, "test", "genome_binner.pkl"), "rb") as r:
#         binners_dict = pickle.load(r)

#     cb = CoolerBinner(cooler_fn, binner=binners_dict[1000000])
#     cooler_binned = cb.to_cooler(outfile)
#     assert os.path.exists(cooler_binned)


# def test_binning_cli_conversion_table():
    
#     infile = os.path.join(dir_test, "test", "cooler.Slc25A37.fragments.hdf5")
#     output_prefix = os.path.join(dir_test, "test", "cli_cooler_binned.hdf5")

#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "reporters",
#             "store",
#             "bins",
#             infile,
#             "-b",
#             "1000000",
#             "-o",
#             output_prefix,
#             "--normalise",
#             "-p",
#             "4",
#             "--conversion_tables",
#             os.path.join(dir_test, "test", "genome_binner.pkl")
#         ],
#     )

#     assert result.exit_code == 0


