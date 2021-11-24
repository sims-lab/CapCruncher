import os
import sys
import subprocess
import shutil
import glob
import pytest



# Pre-run setup
dir_test = os.path.realpath(os.path.dirname(__file__))
dir_test_run = os.path.join(dir_test, 'pipeline_test_run')
dir_package = os.path.dirname(dir_test)
dir_repo = os.path.dirname(dir_package)
dir_pipeline = os.path.join(dir_package, 'pipeline')
dir_data = os.path.join(dir_package, "data")

# data paths
data_path_oligos = os.path.join(dir_data, 'test', 'mm9_capture_oligos_Slc25A37.bed')
data_path_genome = os.path.join(dir_data, 'test', 'data_for_pipeline_run', 'chr14.fa.gz')
data_path_chromsizes = os.path.join(dir_data, 'test', 'data_for_pipeline_run', 'chr14.fa.fai')
data_path_bowtie2_indicies = os.path.join(dir_data, 'test', 'data_for_pipeline_run', 'chr14_bowtie2_indicies')
data_path_config = os.path.join(dir_repo, 'config.yml')


@pytest.fixture(scope="module", autouse=True)
def setup():
    ## Set-up run.
    dir_root = os.getcwd()
    os.mkdir(dir_test_run) 
    os.chdir(dir_test_run)
    cwd = os.getcwd()

    ## Copy fastq file
    for fq in glob.glob(os.path.join(dir_data, 'test', 'data_for_pipeline_run', 'Slc25A37*.fastq.gz')):
        shutil.copy(fq, '.')
    
    # Make indicies
    indicies_dir = os.path.join(dir_data, 'test','data_for_pipeline_run', 'bt2')


    if not os.path.exists(indicies_dir):
        os.mkdir(indicies_dir)
        cmd = f"bowtie2-build {data_path_genome} {indicies_dir}/chr14 --threads 8"
        subprocess.run(cmd.split())
        

    ## Read config and replace with correct paths
    replacements = {'PATH_TO_VIEWPOINTS': data_path_oligos,
                    'PATH_TO_GENOME_FASTA': data_path_genome,
                    'PATH_TO_ALIGNER_INDICIES': f'{indicies_dir}/chr14',
                    'PATH_TO_CHROMOSOME_SIZES': data_path_chromsizes,
                    'PATH_TO_HUB_DIRECTORY': dir_test_run,
                    'PATH_TO_PLOTTING_COORDINATES': os.path.join(dir_data, "test", 'data_for_pipeline_run', 'plot_coords.bed'),
                    'PATH_TO_TSV_FORMATTED_DESIGN_MATRIX': os.path.join(dir_data, "test", 'data_for_pipeline_run', 'design_matrix.tsv'),
                    'PATH_TO_GENES_IN_BED12_FORMAT': os.path.join(dir_data, "test", 'data_for_pipeline_run', 'mm9_chr14_genes.bed'),
                    'HUB_NAME': 'capturec_test',
                    'REGIONS_FOR_NORM': os.path.join(dir_data, "test", 'data_for_pipeline_run', 'regions_for_norm.bed'),
                    }
                
    with open(data_path_config, 'r') as config:
        with open('config.yml', 'w') as writer:
            for line in config:
                for key in replacements:
                    if key in line:
                        line = line.replace(key, replacements[key])

                writer.write(line)
    
    yield

    os.chdir(dir_root)
    shutil.rmtree(dir_test_run)

def test_pipeline_handles_no_drmaa():
    cmd = f'capcruncher pipeline make annotate_sort_viewpoints'
    completed = subprocess.run(cmd.split())
    assert completed.returncode == 0

def test_pipeline_fastq_preprocessing():

    cmd = f'python {dir_pipeline}/pipeline.py make fastq_preprocessing --local -p 4'
    completed = subprocess.run(cmd.split())
    assert completed.returncode == 0
    assert len(glob.glob('capcruncher_preprocessing/digested/*.fastq.gz')) == 12

def test_pipeline_pre_annotation():

    cmd = f'python {dir_pipeline}/pipeline.py make pre_annotation --local -p 4'
    completed = subprocess.run(cmd.split())
    assert completed.returncode == 0

def test_pipeline_post_annotation():
    cmd = f'python {dir_pipeline}/pipeline.py make post_annotation --local -p 4'
    completed = subprocess.run(cmd.split())
    assert completed.returncode == 0

def test_pipeline_post_capcruncher_analysis():
    cmd = f'python {dir_pipeline}/pipeline.py make post_capcruncher_analysis --local -p 4'
    completed = subprocess.run(cmd.split())
    breakpoint()
    assert completed.returncode == 0


def test_pipeline_all():

    cmd = f'python {dir_pipeline}/pipeline.py make full --local -p 4'
    completed = subprocess.run(cmd.split())
    breakpoint()
    assert completed.returncode == 0

def test_stats_exist():

    assert os.path.exists('capcruncher_statistics/capcruncher_statistics.html')

def test_bigwigs_exist():

    assert len(glob.glob('capcruncher_analysis/bigwigs/*.bigWig')) == 16

def test_hub_exists():
    assert os.path.exists('capturec_test.hub.txt')

def test_plot_template_exists():
    try:
        import coolbox
        assert os.path.exists("capcruncher_plots/templates/Slc25A37.pileup.yml")
    except ImportError:
        pass

def test_plot_exists():
    try:
        import coolbox
        assert os.path.exists("capcruncher_plots/Slc25A37_chr14:69878554-69933221.svg")
    except ImportError:
        pass









