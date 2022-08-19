import os
import subprocess
import shutil
import glob
import pytest
import logging


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
        try:
            import requests
            import tarfile

            url = "https://userweb.molbiol.ox.ac.uk/public/asmith/capcruncher/test_indicies.tar.gz"
            output = os.path.join(data_path, "test_indicies.tar.gz")

            r = requests.get(url, stream=True)
            with open(output, 'wb') as f:
                f.write(r.content)

            tar = tarfile.open(output)
            tar.extractall(path=data_path)
            tar.close()
            os.remove(output)
            os.rename(data_path + "/test_indicies", indicies)
            logging.info("Downloaded indicies")

        except Exception as e:
            print(e)
            print("Could not download indicies so generating them")
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
def run_directory_capture(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data_capture")
    return fn

@pytest.fixture(scope="module")
def run_directory_tri(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data_tri")
    return fn

@pytest.fixture(scope="module")
def run_directory_tiled(tmpdir_factory):
    fn = tmpdir_factory.mktemp("data_tiled")
    return fn


@pytest.fixture(scope="module")
def setup_pipeline_run_capture(data_path, run_directory_capture, genome, indicies, config_yaml):

    oligos = os.path.join(data_path, "mm9_capture_oligos_Slc25A37.bed")
    chromsizes = os.path.join(data_path, "chr14.fa.fai")
    fastq = glob.glob(os.path.join(data_path, "*.fastq*"))

    os.chdir(run_directory_capture)

    for fn in fastq:
        shutil.copy(fn, ".")

    ## Read config and replace with correct paths
    replacements = {
        "ANALYSIS_METHOD": "capture",
        "PATH_TO_VIEWPOINTS": oligos,
        "PATH_TO_GENOME_FASTA": genome,
        "PATH_TO_ALIGNER_INDICIES": indicies,
        "PATH_TO_CHROMOSOME_SIZES": chromsizes,
        "PATH_TO_HUB_DIRECTORY": os.path.join(run_directory_capture, "hub_directory"),
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


@pytest.fixture(scope="module")
def setup_pipeline_run_tri(data_path, run_directory_tri, genome, indicies, config_yaml):

    oligos = os.path.join(data_path, "mm9_capture_oligos_Slc25A37.bed")
    chromsizes = os.path.join(data_path, "chr14.fa.fai")
    fastq = glob.glob(os.path.join(data_path, "*.fastq*"))

    os.chdir(run_directory_tri)

    for fn in fastq:
        shutil.copy(fn, ".")

    ## Read config and replace with correct paths
    replacements = {
        "ANALYSIS_METHOD": "tri",
        "PATH_TO_VIEWPOINTS": oligos,
        "PATH_TO_GENOME_FASTA": genome,
        "PATH_TO_ALIGNER_INDICIES": indicies,
        "PATH_TO_CHROMOSOME_SIZES": chromsizes,
        "PATH_TO_HUB_DIRECTORY": os.path.join(run_directory_tri, "hub_directory"),
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

@pytest.fixture(scope="module")
def setup_pipeline_run_tiled(data_path, run_directory_tiled, genome, indicies, config_yaml):

    oligos = os.path.join(data_path, "mm9_capture_oligos_Slc25A37.bed")
    chromsizes = os.path.join(data_path, "chr14.fa.fai")
    fastq = glob.glob(os.path.join(data_path, "*.fastq*"))

    os.chdir(run_directory_tiled)

    for fn in fastq:
        shutil.copy(fn, ".")

    ## Read config and replace with correct paths
    replacements = {
        "ANALYSIS_METHOD": "tiled",
        "PATH_TO_VIEWPOINTS": oligos,
        "PATH_TO_GENOME_FASTA": genome,
        "PATH_TO_ALIGNER_INDICIES": indicies,
        "PATH_TO_CHROMOSOME_SIZES": chromsizes,
        "PATH_TO_HUB_DIRECTORY": os.path.join(run_directory_tiled, "hub_directory"),
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

@pytest.mark.order(1)
def test_pipeline_capture(setup_pipeline_run_capture):

    cmd = f"capcruncher pipeline make full --local -p 8"
    completed = subprocess.run(cmd.split())
    assert completed.returncode == 0

@pytest.mark.order(1)
def test_pipeline_tri(setup_pipeline_run_tri):

    cmd = f"capcruncher pipeline make full --local -p 8"
    completed = subprocess.run(cmd.split())
    assert completed.returncode == 0

@pytest.mark.order(1)
def test_pipeline_tiled(setup_pipeline_run_tiled):

    cmd = f"capcruncher pipeline make full --local -p 8"
    completed = subprocess.run(cmd.split())
    assert completed.returncode == 0

@pytest.mark.order(2)
def test_digested_exists(run_directory_capture):
    assert len(
        glob.glob(f"{run_directory_capture}/capcruncher_preprocessing/digested/*.fastq*")
    ) == (4 * 2)

@pytest.mark.order(2)
def test_stats_exist(run_directory_capture):
    assert os.path.exists(
        f"{run_directory_capture}/capcruncher_statistics/capcruncher_statistics.html"
    )

@pytest.mark.order(2)
@pytest.mark.parametrize("n_samples,n_groups,n_viewpoints", [(4, 2, 1)])
def test_bigwigs_exist(run_directory_capture, n_samples, n_groups, n_viewpoints):
    import math

    n_bigwigs_expected = sum(
        [
            (n_samples * len(["raw", "normalised"]) * n_viewpoints),
            (n_groups * n_viewpoints),
            (math.perm(n_groups, 2) * n_viewpoints),
        ],
    )
    assert (
        len(glob.glob(f"{run_directory_capture}/capcruncher_analysis/bigwigs/*.bigWig"))
        == n_bigwigs_expected
    )

@pytest.mark.order(2)
def test_hub_exists(run_directory_capture):
    assert os.path.exists(f"{run_directory_capture}/hub_directory")

@pytest.mark.order(2)
def test_plot_template_exists(run_directory_capture):
    try:
        import coolbox

        assert os.path.exists(
            f"{run_directory_capture}/capcruncher_plots/templates/Slc25A37.pileup.yml"
        )
    except ImportError:
        pass

@pytest.mark.order(2)
def test_plot_exists(run_directory_capture):
    try:
        import coolbox

        assert os.path.exists(
            f"{run_directory_capture}/capcruncher_plots/Slc25A37_chr14:69878554-69933221.svg"
        )
    except ImportError:
        pass
