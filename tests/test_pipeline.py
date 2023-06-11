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


@pytest.fixture(scope="module")
def temp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("test_pipeline")


@pytest.fixture(scope="module")
def package_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    return os.path.dirname(dirname)


@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "data_for_pipeline_run")
    return data_dir


@pytest.fixture(scope="module")
def fasta(data_path):
    return os.path.join(data_path, "chr14.fa.gz")


@pytest.fixture(scope="module")
def genome():
    return "mm9"


@pytest.fixture(scope="module")
def indicies(data_path, genome):
    indicies = os.path.join(data_path, "chr14_bowtie2_indicies")
    if not os.path.exists(indicies):
        try:
            import requests
            import tarfile

            url = "https://userweb.molbiol.ox.ac.uk/public/asmith/capcruncher/test_indicies.tar.gz"
            output = os.path.join(data_path, "test_indicies.tar.gz")

            r = requests.get(url, stream=True)
            with open(output, "wb") as f:
                f.write(r.content)

            tar = tarfile.open(output)
            tar.extractall(path=data_path)
            tar.close()
            os.remove(output)
            os.rename(data_path + "/test_indicies", indicies)
            logger.info("Downloaded indicies")

        except Exception as e:
            print(e)
            print("Could not download indicies so generating them")
            os.mkdir(indicies)
            cmd = f"bowtie2-build {genome} {indicies}/bt2 --threads 8"
            subprocess.run(cmd.split())

    return os.path.join(indicies, "bt2")


@pytest.fixture(scope="module")
def binsizes():
    return np.random.randint(int(1e3), int(1e6), size=3)


@pytest.fixture(scope="module")
def plot_coords(data_path):
    return os.path.join(data_path, "mm9_capture_viewpoints_Slc25A37.bed")


@pytest.fixture(scope="module")
def design(data_path):
    return os.path.join(data_path, "design_matrix.tsv")


@pytest.fixture(scope="module")
def viewpoints(data_path):
    return os.path.join(data_path, "mm9_capture_viewpoints_Slc25A37.bed")


@pytest.fixture(scope="module")
def fastqs(data_path):
    return glob.glob(os.path.join(data_path, "*.fastq*"))


@pytest.fixture(scope="module")
def chromsizes(data_path):
    return os.path.join(data_path, "chr14.fa.fai")


@pytest.fixture(scope="module")
def run_dir_capture(temp_dir):
    current_date = datetime.now().strftime("%Y-%m-%d")
    project_id = "project_name"
    assay = "capture"
    return os.path.join(temp_dir, f"{current_date}_{project_id}_{assay}")


@pytest.fixture(scope="module")
def run_dir_tiled(temp_dir):
    current_date = datetime.now().strftime("%Y-%m-%d")
    project_id = "project_name"
    assay = "tiled"
    return os.path.join(temp_dir, f"{current_date}_{project_id}_{assay}")


@pytest.fixture(scope="module", params=["capture", "tri", "tiled"])
def config(
    temp_dir,
    package_path,
    fasta,
    indicies,
    binsizes,
    viewpoints,
    chromsizes,
    design,
    plot_coords,
    request,
):
    cwd = os.getcwd()
    os.chdir(temp_dir)

    method = request.param

    cookiecutter(
        f"{package_path}/pipeline/config/",
        extra_context={
            "method": method,
            "design": design,
            "viewpoints": viewpoints,
            "genome": "mm9",
            "is_custom_genome": "no",
            "genome_organism": "Mus musculus",
            "genome_fasta": fasta,
            "genome_chromosome_sizes": chromsizes,
            "genome_indicies": indicies,
            "restriction_enzyme": "dpnii",
            "remove_blacklist": "no",
            "genomic_bin_size": binsizes,
            "prioritize_cis_slices": "yes",
            "priority_chromosomes": "viewpoints",
            "make_ucsc_hub": "yes",
            "ucsc_hub_directory": ".",
            "ucsc_hub_name": "CCHUB_TEST",
            "ucsc_hub_email": "Email address (UCSC required)",
            "ucsc_track_color_by": "samplename",
            "make_plots": "yes",
            "plotting_coordinates": plot_coords,
            "plotting_normalisation": "n_interactions",
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
def test_pipeline(config):
    from sh import capcruncher

    capcruncher.pipeline("-c", "8", "all", "-p")


@pytest.mark.order(2)
def test_stats_exist(run_dir_capture):
    assert os.path.exists(
        f"{run_dir_capture}/capcruncher_output/statistics/capcruncher_report.html"
    )


@pytest.mark.order(2)
@pytest.mark.parametrize("n_samples,n_groups,n_viewpoints", [(4, 2, 1)])
def test_bigwigs_exist(run_dir_capture, n_samples, n_groups, n_viewpoints):
    import math

    n_bigwigs_expected = sum(
        [
            (n_samples * len(["raw", "normalised"]) * n_viewpoints),
            (n_groups * n_viewpoints),
            (math.perm(n_groups, 2) * n_viewpoints),
        ],
    )

    bigwigs = list(
        pathlib.Path(f"{run_dir_capture}/capcruncher_output/pileups/bigwigs/").glob(
            "*.bigWig"
        )
    )
    assert len(bigwigs) == n_bigwigs_expected


@pytest.mark.order(2)
def test_reporters_are_binned(run_dir_tiled, binsizes):
    import cooler

    example_cooler = os.path.join(
        run_dir_tiled, "capcruncher_output/pileups/counts/SAMPLE-A_REP1.hdf5"
    )
    cooler_groups = cooler.api.list_coolers(example_cooler)
    assert len(cooler_groups) == len(binsizes) + 1


@pytest.mark.order(2)
def test_hub_exists(run_dir_capture):
    assert os.path.exists(f"{run_dir_capture}/CCHUB_TEST")
