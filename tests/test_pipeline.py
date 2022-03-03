import os
import subprocess
import shutil
import glob
import pytest


@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "data_for_pipeline_run")
    return data_dir


@pytest.fixture(scope="module")
def genome(data_path):
    return os.path.join(data_path, "chr14.fa.gz")


@pytest.fixture(scope="module")
def indicies(data_path, genome):

    indicies = os.path.join(data_path, "chr14_bowtie2_indicies")
    if not os.path.exists(indicies):
        os.mkdir(indicies)
        cmd = f"bowtie2-build {genome} {indicies}/bt2 --threads 8"
        subprocess.run(cmd.split())

    return os.path.join(indicies, "bt2")


@pytest.fixture(scope="module")
def config_yaml(data_path):
    current_dir = os.path.dirname(os.path.realpath(__file__))
    repo_dir = os.path.dirname(current_dir)
    config = os.path.join(repo_dir, "config.yml")
    return config


@pytest.fixture(scope="module")
def run_directory(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data")
    return fn


@pytest.fixture(scope="module")
def setup_pipeline_run(data_path, run_directory, genome, indicies, config_yaml):

    oligos = os.path.join(data_path, "mm9_capture_oligos_Slc25A37.bed")
    chromsizes = os.path.join(data_path, "chr14.fa.fai")
    fastq = glob.glob(os.path.join(data_path, "*.fastq*"))

    os.chdir(run_directory)

    for fn in fastq:
        shutil.copy(fn, ".")

    ## Read config and replace with correct paths
    replacements = {
        "PATH_TO_VIEWPOINTS": oligos,
        "PATH_TO_GENOME_FASTA": genome,
        "PATH_TO_ALIGNER_INDICIES": indicies,
        "PATH_TO_CHROMOSOME_SIZES": chromsizes,
        "PATH_TO_HUB_DIRECTORY": os.path.join(run_directory, "hub_directory"),
        "PATH_TO_PLOTTING_COORDINATES": os.path.join(data_path, "plot_coords.bed"),
        "PATH_TO_TSV_FORMATTED_DESIGN_MATRIX": os.path.join(
            data_path, "design_matrix.tsv"
        ),
        "PATH_TO_GENES_IN_BED12_FORMAT": os.path.join(data_path, "mm9_chr14_genes.bed"),
        "HUB_NAME": "CAPCRUNCHER_TEST_HUB",
        "REGIONS_FOR_NORM": os.path.join(data_path, "regions_for_norm.bed"),
    }

    with open(config_yaml, "r") as config:
        with open("config.yml", "w") as writer:
            for line in config:
                for key in replacements:
                    if key in line:
                        line = line.replace(key, replacements[key])

                writer.write(line)

    yield


# def test_pipeline_handles_no_drmaa(setup_pipeline_run):
#     cmd = f'capcruncher pipeline make annotate_sort_viewpoints'
#     completed = subprocess.run(cmd.split())
#     assert completed.returncode == 0


def test_pipeline_all(setup_pipeline_run):

    cmd = f"capcruncher pipeline make full --local -p 8"
    completed = subprocess.run(cmd.split())
    assert completed.returncode == 0


def test_digested_exists(run_directory):
    assert len(
        glob.glob(f"{run_directory}/capcruncher_preprocessing/digested/*.fastq*")
    ) == (4 * 2)


def test_stats_exist(run_directory):
    assert os.path.exists(
        f"{run_directory}/capcruncher_statistics/capcruncher_statistics.html"
    )


@pytest.mark.parametrize("n_samples,n_groups,n_viewpoints", [(4, 2, 1)])
def test_bigwigs_exist(run_directory, n_samples, n_groups, n_viewpoints):
    import math

    n_bigwigs_expected = sum(
        [
            (n_samples * len(["raw", "normalised"]) * n_viewpoints),
            (n_groups * n_viewpoints),
            (math.perm(n_groups, 2) * n_viewpoints),
        ],
    )
    assert (
        len(glob.glob(f"{run_directory}/capcruncher_analysis/bigwigs/*.bigWig"))
        == n_bigwigs_expected
    )


def test_hub_exists(run_directory):
    assert os.path.exists(f"{run_directory}/hub_directory")


def test_plot_template_exists(run_directory):
    try:
        import coolbox

        assert os.path.exists(
            f"{run_directory}/capcruncher_plots/templates/Slc25A37.pileup.yml"
        )
    except ImportError:
        pass


def test_plot_exists(run_directory):
    try:
        import coolbox

        assert os.path.exists(
            f"{run_directory}/capcruncher_plots/Slc25A37_chr14:69878554-69933221.svg"
        )
    except ImportError:
        pass
