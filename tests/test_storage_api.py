import pyranges1 as pr
import pandas as pd

from capcruncher.api.storage import CoolerBinner, Viewpoint


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
