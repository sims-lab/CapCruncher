# import os
# import sys
# import click
# from pybedtools.bedtool import BedTool

# from capcruncher.tools.pileup import CoolerBedGraph

# # Pre-run setup
# dir_test = os.path.realpath(os.path.dirname(__file__))
# dir_package = os.path.dirname(dir_test)
# dir_data = os.path.join(dir_package, "data")


# def test_raw_bedgraph():

#     clr = os.path.join(dir_data, "test", "test.Slc25A37.hdf5")
#     ccbdg = CoolerBedGraph(clr)
#     ccbdg.normalise_bedgraph(method="raw")

# def test_n_cis_bedgraph():

#     clr = os.path.join(dir_data, "test", "test.Slc25A37.hdf5")
#     ccbdg = CoolerBedGraph(clr)
#     ccbdg.normalise_bedgraph(method="n_cis")

# def test_region_bedgraph():

#     clr = os.path.join(dir_data, "test", "test.Slc25A37.hdf5")
#     ccbdg = CoolerBedGraph(clr)

#     viewpoints = BedTool(os.path.join(dir_data, "test", "mm9_capture_oligos_Slc25A37.bed"))
#     region = viewpoints.slop(b=2e5, genome="mm9")
#     region = region.subtract(viewpoints)
#     ccbdg.normalise_bedgraph(method="region", region=region.fn)



    

