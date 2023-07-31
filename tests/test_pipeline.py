import os
import subprocess
import shutil
import glob
import pytest
from loguru import logger
import numpy as np
import pathlib
from cookiecutter.main import cookiecutter
from datetime import datetime


# Fixtures
@pytest.fixture(scope="module")
def test_dir(tmpdir_factory):
    return pathlib.Path(tmpdir_factory.mktemp("test_dir"))


@pytest.fixture(scope="module")
def repo_path():
    fn = pathlib.Path(__file__).resolve()
    dirname = fn.parent
    return dirname.parent


@pytest.fixture(scope="module")
def package_path(repo_path):
    return repo_path.joinpath("capcruncher")


@pytest.fixture(scope="module")
def data_path():
    fn = pathlib.Path(__file__).resolve()
    dirname = fn.parent
    data_dir = dirname.joinpath("data", "data_for_pipeline_run")
    return data_dir


@pytest.fixture(scope="module")
def fasta(data_path):
    return data_path.joinpath("chr14.fa.gz")


@pytest.fixture(scope="module")
def indicies(data_path, genome):
    indicies = data_path.joinpath("chr14_bowtie2_indicies")
    if not indicies.exists():
        try:
            import requests
            import tarfile

            url = "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/capcruncher/test_indicies.tar.gz"
            output = data_path.joinpath("test_indicies.tar.gz")

            r = requests.get(url, stream=True)
            with output.open("wb") as f:
                f.write(r.content)

            tar = tarfile.open(output)
            tar.extractall(path=data_path)
            tar.close()
            output.unlink()
            (data_path / "test_indicies").rename(indicies)
            logger.info("Downloaded indicies")

        except Exception as e:
            print(e)
            print("Could not download indicies so generating them")
            indicies.mkdir()
            cmd = f"bowtie2-build {genome} {indicies}/bt2 --threads 8"
            subprocess.run(cmd.split())

    return indicies.joinpath("bt2")


@pytest.fixture(scope="module")
def plot_coords(data_path):
    return data_path.joinpath("plot_coords.bed")


@pytest.fixture(scope="module")
def design(data_path):
    return data_path.joinpath("design_matrix.tsv")


@pytest.fixture(scope="module")
def viewpoints(data_path):
    return data_path.joinpath("mm9_capture_viewpoints_Slc25A37.bed")


@pytest.fixture(scope="module")
def fastqs(data_path):
    return list(data_path.glob("*.fastq*"))


@pytest.fixture(scope="module")
def chromsizes(data_path):
    return data_path.joinpath("chr14.fa.fai")


@pytest.fixture(scope="module")
def run_dir_capture(test_dir):
    current_date = datetime.now().strftime("%Y-%m-%d")
    project_id = "project_name"
    assay = "capture"
    return test_dir.joinpath(f"{current_date}_{project_id}_{assay}")


@pytest.fixture(scope="module")
def run_dir_tiled(test_dir):
    current_date = datetime.now().strftime("%Y-%m-%d")
    project_id = "project_name"
    assay = "tiled"
    return test_dir.joinpath(f"{current_date}_{project_id}_{assay}")


@pytest.fixture(scope="module")
def genome():
    return "mm9"


@pytest.fixture(scope="module")
def binsizes():
    return [10000, 20000, 40000, 80000, 160000, 320000, 640000]


@pytest.fixture(scope="module")
def hub_dir(run_dir_capture):
    return run_dir_capture / "HUB_DIR"


@pytest.fixture(scope="module", params=["capture", "tri", "tiled"])
def config(
    test_dir,
    package_path,
    fasta,
    fastqs,
    indicies,
    binsizes,
    viewpoints,
    chromsizes,
    design,
    plot_coords,
    request,
):
    cwd = pathlib.Path.cwd()
    os.chdir(test_dir)

    METHODS = {"capture": "Capture-C", "tri": "Tri-C", "tiled": "Tiled-C"}
    method = METHODS[request.param]

    cookiecutter(
        f"{package_path}/pipeline/config/",
        extra_context={
            "method": str(method),
            "design": str(design),
            "viewpoints": str(viewpoints),
            "genome": "mm9",
            "is_custom_genome": "no",
            "genome_organism": "Mus musculus",
            "genome_fasta": str(fasta),
            "genome_chromosome_sizes": str(chromsizes),
            "genome_indicies": str(indicies),
            "restriction_enzyme": "dpnii",
            "remove_blacklist": "no",
            "genomic_bin_size": " ".join([str(b) for b in binsizes]),
            "prioritize_cis_slices": "yes",
            "priority_chromosomes": "viewpoints",
            "make_ucsc_hub": "yes",
            "ucsc_hub_directory": "HUB_DIR",
            "ucsc_hub_name": "CCHUB_TEST",
            "ucsc_hub_email": "Email address (UCSC required)",
            "ucsc_track_color_by": "samplename",
            "make_plots": "yes",
            "plotting_coordinates": str(plot_coords),
            "plotting_normalisation": "n_interactions",
            "differential_contrast": "condition",
        },
        no_input=True,
    )

    # Move config files and fastq files
    current_date = datetime.now().strftime("%Y-%m-%d")
    project_id = "project_name"
    assay = method.lower().split("-")[0]
    os.chdir(f"{current_date}_{project_id}_{assay}")

    for fq in fastqs:
        fq_new = pathlib.Path(fq.name)
        fq_new.symlink_to(fq)

    yield

    os.chdir(cwd)


@pytest.mark.order(1)
def test_pipeline(config, cores):
    from sh import capcruncher
    import sys
    import subprocess

    if cores:
        cores = cores
    else:
        cores = 1

    try:
        result = subprocess.run(
            ["capcruncher", "pipeline", "-c", str(cores), "all", "-p"]
        )
    except Exception as e:
        print(e)
        raise e

    assert result.returncode == 0


@pytest.mark.order(2)
def test_stats_exist(run_dir_capture):
    run_dir_capture = pathlib.Path(run_dir_capture)
    assert (
        run_dir_capture / "capcruncher_output/results/capcruncher_report.html"
    ).exists()


@pytest.mark.order(2)
@pytest.mark.parametrize("n_samples,n_groups,n_viewpoints", [(4, 2, 1)])
def test_bigwigs_exist(run_dir_capture, n_samples, n_groups, n_viewpoints):
    import math

    run_dir_capture = pathlib.Path(run_dir_capture)
    n_bigwigs_expected = sum(
        [
            (n_samples * len(["raw", "normalised"]) * n_viewpoints),
            (n_groups * n_viewpoints),
            (math.perm(n_groups, 2) * n_viewpoints),
        ],
    )

    bigwigs = list(
        pathlib.Path(run_dir_capture / "capcruncher_output/results/").glob(
            "**/*.bigWig"
        )
    )
    assert len(bigwigs) == n_bigwigs_expected


@pytest.mark.order(2)
def test_reporters_are_binned(run_dir_tiled, binsizes):
    import cooler

    example_cooler = (
        run_dir_tiled / "capcruncher_output/results/SAMPLE-A_REP1/SAMPLE-A_REP1.hdf5"
    )
    cooler_groups = cooler.api.list_coolers(str(example_cooler))
    assert len(cooler_groups) == len(binsizes) + 1


@pytest.mark.order(2)
def test_hub_exists(hub_dir):
    assert hub_dir.exists()
