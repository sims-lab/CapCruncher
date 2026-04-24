import pyranges1 as pr

from capcruncher.api.storage import Viewpoint


def test_viewpoint_from_bed_returns_pyranges(tmp_path):
    bed_path = tmp_path / "viewpoints.bed"
    bed_path.write_text("chr1\t10\t20\tcapture_Slc25A37\n")

    viewpoint = Viewpoint.from_bed(
        bed=str(bed_path), viewpoint="Slc25A37", assay="capture"
    )

    assert isinstance(viewpoint.coordinates, pr.PyRanges)
    assert viewpoint.chromosomes == ["chr1"]
