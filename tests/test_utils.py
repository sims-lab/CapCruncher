import os
import itertools
import pytest
import xopen
import gzip
import subprocess
from click.testing import CliRunner
from capcruncher.cli import cli


@pytest.fixture(scope="module")
def cli_runner():
    return CliRunner()


@pytest.fixture(scope="module")
def data_path_pipeline():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "data_for_pipeline_run")
    return data_dir

@pytest.fixture(scope="module")
def data_path_utils():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "utils")
    return data_dir


@pytest.fixture(scope="module")
def genome(data_path_pipeline):
    return os.path.join(data_path_pipeline, "chr14.fa.gz")


@pytest.fixture(scope="module")
def indicies(data_path_pipeline, genome):

    indicies = os.path.join(data_path_pipeline, "chr14_bowtie2_indicies")
    if not os.path.exists(indicies):
        os.mkdir(indicies)
        cmd = f"bowtie2-build {genome} {indicies}/bt2 --threads 8"
        subprocess.run(cmd.split())

    return os.path.join(indicies, "bt2")


@pytest.mark.parametrize(
    "viewpoints,flags",
    [
        (
            "viewpoints.fa",
            ["-r", "dpnii"],
        ),
    ],
)
def test_viewpoint_coordinates(cli_runner, viewpoints, genome, indicies, data_path_pipeline, data_path_utils, flags, tmpdir):

    viewpoints = os.path.join(data_path_utils, viewpoints)
    outfile = os.path.join(tmpdir, "viewpoint_coords.bed")

    result = cli_runner.invoke(
        cli, ["utilities", 
              "viewpoint-coordinates", 
              "-v", 
              viewpoints,
              "-i",
              indicies,
              "-g",
              genome,
              "-o", 
              outfile,  
              *flags]
    )
    assert result.exit_code == 0
    assert os.path.exists(outfile)
