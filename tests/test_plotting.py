import pathlib

from plotnado import GenomicFigure


def test_plotnado_renders_representative_tracks(tmp_path):
    data_path = pathlib.Path(__file__).resolve().parent / "data"
    bigwig = (
        data_path / "test_bigwigs" / "Slc25A37-test-1x_1.normalised.Slc25A37.bigWig"
    )
    bed = data_path / "data_for_pipeline_run" / "mm9_capture_viewpoints_Slc25A37.bed"
    output = tmp_path / "plotnado_tracks.png"

    fig = GenomicFigure()
    fig.scalebar()
    fig.bigwig(str(bigwig), title=bigwig.stem, min_value=0)
    fig.bed(str(bed), title="capture viewpoints")
    fig.axis()
    fig.save(output, region="chr14:69868303-69956880")

    assert output.exists()


def test_plotnado_toml_round_trip(tmp_path):
    data_path = pathlib.Path(__file__).resolve().parent / "data"
    bigwig = (
        data_path / "test_bigwigs" / "Slc25A37-test-1x_1.normalised.Slc25A37.bigWig"
    )
    toml_path = tmp_path / "plotnado_template.toml"

    fig = GenomicFigure()
    fig.bigwig(str(bigwig), title="test")
    fig.add_track(
        "capcruncher",
        file=str(data_path / "reporters_store" / "SAMPLE-A_REP1_binned.hdf5"),
        title="Slc25A37",
        resolution=5000,
        viewpoint="Slc25A37",
        normalisation="raw",
        balance=False,
    )
    fig.to_toml(toml_path)

    fig2 = GenomicFigure.from_toml(toml_path)

    assert toml_path.exists()
    assert [track.__class__.__name__ for track in getattr(fig2, "tracks")] == [
        "BigWigTrack",
        "CapcruncherTrack",
    ]
