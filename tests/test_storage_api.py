import pyranges1 as pr
import pandas as pd
import pytest

from capcruncher.api.interactions.cooler.binning import CoolerBinner
from capcruncher.api.interactions.cooler.create import create_cooler_cc
from capcruncher.api.interactions.cooler.viewpoints import Viewpoint


def test_create_cooler_rejects_empty_pixel_table(tmp_path):
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [0],
            "end": [100],
            "name": [0],
        }
    )
    viewpoints = tmp_path / "viewpoints.bed"
    viewpoints.write_text("chr1\t0\t100\tvp1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="empty pixels table"):
        create_cooler_cc(
            tmp_path / "empty_counts.hdf5",
            bins=bins,
            pixels=pd.DataFrame(),
            viewpoint_name="vp1",
            viewpoint_path=viewpoints,
        )


def test_viewpoint_from_bed_returns_pyranges(tmp_path):
    bed_path = tmp_path / "viewpoints.bed"
    bed_path.write_text("chr1\t10\t20\tcapture_Slc25A37\n")

    viewpoint = Viewpoint.from_bed(
        bed=str(bed_path), viewpoint="Slc25A37", assay="capture"
    )

    assert isinstance(viewpoint.coordinates, pr.PyRanges)
    assert viewpoint.chromosomes == ["chr1"]


def test_viewpoint_from_bed_matches_literal_viewpoint_name(tmp_path):
    bed_path = tmp_path / "viewpoints.bed"
    bed_path.write_text(
        "chr1\t10\t20\tcapture_Slc25A37\n"
        "chr1\t30\t40\tcapture_Slc25A37_extra\n"
    )

    viewpoint = Viewpoint.from_bed(
        bed=str(bed_path), viewpoint="Slc25A37", assay="capture"
    )

    assert viewpoint.coords == ["chr1:10-20"]


def test_midpoint_fragment_mapping_uses_integer_coordinates():
    binner = CoolerBinner.__new__(CoolerBinner)
    binner.method = "midpoint"
    binner.minimum_overlap = 0.51
    binner.__dict__["genomic_bins"] = pr.PyRanges(
        pd.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [0],
                "End": [10],
                "genomic_bin_id": [0],
            }
        )
    )
    binner.__dict__["fragment_bins"] = pr.PyRanges(
        pd.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [0],
                "End": [9],
                "fragment_id": [1],
            }
        )
    )

    mapped = binner.fragment_to_genomic_table

    assert mapped["Start_b"].tolist() == [4]
    assert mapped["End_b"].tolist() == [5]
