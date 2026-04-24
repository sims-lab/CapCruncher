# Pipeline

The CapCruncher pipeline handles the processing of raw data from the sequencer to the generation of a contact matrix, generation of plots and production of a UCSC genome browser track hub.

This pipeline is based on the Snakemake workflow management system. Snakemake is a Python-based workflow management system that allows for the creation of reproducible and scalable data analyses. All elements of the workflow have been wrapped into the CapCruncher Python package. This allows for the pipeline to be run using the `capcruncher pipeline` command rather than having to run the pipeline using Snakemake directly.

Checkout the [Hints and Tips](tips.md) page for some useful tips on configuring and running the pipeline.

## Pipeline Configuration

### Configuration File

The pipeline is configured using a YAML file. It is strongly recommended to use the `capcruncher pipeline-config` command to generate a template configuration file. This command will generate a template configuration file with all available options and descriptions of each option.

``` bash
capcruncher pipeline-config
```

This utility will walk through the configuration options and generate a configuration file. It will generate a new directory <DATE>_<PROJECT NAME>_<ASSAY> and place the filled-out `capcruncher_config.yml` file in this directory.

For an example configuration file, see [here](examples/capcruncher_config.yml).

The configuration file can be edited manually if required e.g. to add a manually generated `design` file. Just ensure that the configuration file is valid YAML. A common error is to use tabs instead of spaces, this will cause the pipeline to fail while parsing the configuration file.

All options in the configuration file are documented within the file itself. Only the required options need to be filled out. The pipeline will use default values for all other options.

### Design File

The design file is a tab/comma/space-delimited file that contains the sample names and the metadata for each sample. This file is completely optional and only used for comparisons between Capture-C and Tri-C data. If it is not provided the pipeline will perform a basic sample name comparison to generate a basic design file. However, this will not be as accurate as a manually generated design file. The `design` file is a tab delimited file with the following columns:

- `sample`: The name of the FASTQ file (without the _R1.fastq.gz or_2.fastq.gz suffix)
- `condition`: The Group that the sample belongs to.

Provide the path to this file in the config file under the `design` key.

### Setting up the input directory

Ensure that you have the [configuration file](#configuration-file) and the fastq files to process in the current working directory. Symbolic links can be used to link to the fastq files if they are stored elsewhere but please ensure that the full path to the fastq files is used to create the symbolic links. e.g.

``` bash
# Example 1 - Make a symbolic link to the fastq file in the current directory
ln -s /ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/SAMPLE-A_REP1_1.fastq.gz .

# Example 2 - Make a symbolic link to the fastq file in the current directory with a different name
ln -s /ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/SAMPLE-A_REP1_1.fastq.gz SAMPLE-A_REP1_1.fastq.gz

# Example 3 - Use realpath to get the full path to the fastq file and then create a symbolic link to it in another directory
ln -s $(realpath SAMPLE-A_REP1_1.fastq.gz) /tmp/pytest-of-asmith/pytest-current/test_dircurrent/
```

The pipeline will automatically detect the configuration file and the fastq files. For example your working directory should look like this:

``` bash
2023-08-02_project_name_capture/
|-- SAMPLE-A_REP1_1.fastq.gz -> /ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/SAMPLE-A_REP1_1.fastq.gz
|-- SAMPLE-A_REP1_2.fastq.gz -> /ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/SAMPLE-A_REP1_2.fastq.gz
|-- SAMPLE-A_REP2_1.fastq.gz -> /ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/SAMPLE-A_REP2_1.fastq.gz
|-- SAMPLE-A_REP2_2.fastq.gz -> /ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/SAMPLE-A_REP2_2.fastq.gz
|-- SAMPLE-B_REP1_1.fastq.gz -> /ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/SAMPLE-B_REP1_1.fastq.gz
|-- SAMPLE-B_REP1_2.fastq.gz -> /ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/SAMPLE-B_REP1_2.fastq.gz
|-- SAMPLE-B_REP2_1.fastq.gz -> /ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/SAMPLE-B_REP2_1.fastq.gz
|-- SAMPLE-B_REP2_2.fastq.gz -> /ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/SAMPLE-B_REP2_2.fastq.gz
`-- capcruncher_config.yml
```

## Running the pipeline

### Basic Usage

The pipeline is run using the `capcruncher pipeline` command.

``` bash
# Usage
capcruncher pipeline --cores <NUMBER OF CORES TO USE>

# Example
capcruncher pipeline --cores 8
```

### HPC Cluster Usage (Recommended if available)

The pipeline can also be run on HPC clusters using a number of different job schedulers. See [here](cluster_config.md) for a quick overview on how to configure the pipeline to run on HPC clusters.

For further information see both the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) and the [Snakemake profile documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).

This is a quick example of how to run the pipeline with a pre-generated profile. This is not a complete guide and you will need to modify the configuration to suit your cluster.

``` bash
capcruncher pipeline -c <NUMBER OF CORES e.g. 20> --preset <PRESET NAME>
```

Install bundled editable profiles with `capcruncher pipeline-init`. They are
written to `${XDG_CONFIG_HOME:-~/.config}/snakemake`, where they can be edited
like normal Snakemake profiles. See the [cluster setup guide](cluster_config.md)
for update and customization instructions.

### Apptainer Usage (Recommended if available)

Containers have the advantage of their contents being fixed at the time of creation. This means that the pipeline will always run with the same versions of the software and aids reliablity and reproducibility. The pipeline can be run using Apptainer containers through the bundled Snakemake presets. This is the recommended method of running the pipeline on clusters that provide Apptainer.

Apptainer is the supported container runtime on HPC systems. Docker is intended for local workstation and CI usage, not shared cluster execution.

The pipeline can be run using the bundled Snakemake 9 Apptainer presets.

``` bash
# Local mode
capcruncher pipeline --preset capcruncher-local-apptainer --cores <NUMBER OF CORES TO USE>

# Cluster mode
capcruncher pipeline --preset capcruncher-slurm-apptainer --cores <NUMBER OF CORES TO USE>
```

You can also run the CapCruncher container directly with Apptainer:

``` bash
apptainer exec \
  --bind "$PWD":/work \
  --pwd /work \
  docker://ghcr.io/sims-lab/capcruncher:latest \
  capcruncher pipeline --cores 8
```

### Docker Usage

Docker is supported for running the whole CapCruncher CLI image. This is most useful on local workstations and CI runners. For HPC systems, use Apptainer.

``` bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -v "$PWD":/work \
  -w /work \
  ghcr.io/sims-lab/capcruncher:latest \
  pipeline --cores 8
```

The command must be run from the directory containing `capcruncher_config.yml` and the FASTQ files or mounted symlinks. See the [Docker guide](docker.md) for details.

The Docker image also includes Apptainer for nested containerised workflow execution. If your Docker host permits nested Apptainer, run the image with the required host privileges and use the `capcruncher-local-apptainer` preset inside the container. See the [Docker guide](docker.md) for the full command.

For Snakemake-managed job containers on the host, use the Apptainer presets instead of a Docker preset. Snakemake uses Apptainer to execute `docker://` container image URIs, including the default CapCruncher image configured in `capcruncher_config.yml`.

### Avoiding Disconnection from the Cluster

In order to avoid disconnecting from the cluster, it is recommended to run the pipeline in a [tmux](https://linuxize.com/post/getting-started-with-tmux/) session. Alternatively, [nohup](https://linuxize.com/post/linux-nohup-command/) can be used to run the pipeline in the background. For example:

``` bash
# tmux example
tmux new -s capcruncher
capcruncher pipeline --cores 8 --preset capcruncher-slurm-apptainer

# nohup example
nohup capcruncher pipeline --cores 8 --preset capcruncher-slurm-apptainer &
```

## Pipeline Steps

### All Assays

For all assays the pipeline consists of the following steps:

1. **Quality Control**: FastQC is used to perform quality control on the FASTQ files.
1. **Read Splitting**: The FASTQ files are split into parts of a user-defined size (default 1 million reads per part). This is done to allow for parallel processing of the FASTQ files and to reduce the memory requirements of the pipeline.
1. **Remove PCR Duplicates**: PCR duplicates are removed from the FASTQ files using the CapCruncher package to reduce the memory and CPU requirements of the pipeline.
1. **Read Trimming**: Trimming of the FASTQ files is performed using Trim Galore.
1. **Read Combining**: The trimmed FASTQ files are combined using FLASh to obtain the ligtion junctions.
1. **Read _in silico_ Digestion**: The combined FASTQ files are digested _in silico_ using the restriction enzyme or site specified in the configuration file.
1. **Read Alignment**: The digested FASTQ files are aligned to the reference genome using bowtie2.
1. **Alignment Annotation**: The aligned reads are annotated using the CapCruncher package.
1. **Alignment Filtering**: The aligned reads are filtered using the CapCruncher package.
1. **Alignment PCR Duplicate Removal**: PCR duplicates are removed from the aligned reads using the CapCruncher package.
1. **Contact Matrix Generation**: Contact matrices are generated using the CapCruncher package and stored in cooler (HDF5) format.
1. **Pipeline Statistics**: Statistics are generated for each sample using the CapCruncher package.
1. **Pipeline Plots**: Plots and PlotNado-compatible templates (TOML format) are generated for each sample using the CapCruncher package.

### Capture-C and Tri-C

1. **Pileup Generation**: BigWig files are generated for each sample using the CapCruncher package.
1. **Pileup Normalisation**: The pileup files are normalised using the CapCruncher package.
1. **Pileup Comparison**: Interactions are compared between samples (if two or more replicates are provided) using the CapCruncher package.

!!! note
    Additional methods for comparing interactions between samples will be added in the future.

1. **Differenital Interaction Analysis**: Differential interactions are identified between groups of samples using [PyDESeq2](https://github.com/owkin/PyDESeq2).

!!! warning
    The UCSC Genome Browser Track Hub is only generated if the `ucsc_genome_browser_track_hub` option is set to `True` in the configuration file.

1. **UCSC Genome Browser Track Hub**: A UCSC Genome Browser track hub is generated using the CapCruncher package.

### Tiled-C

1. **Contact Matrix Normalisation**: The contact matrices are normalised using the CapCruncher package with various third-party integrations.
1. **Plot Generation**: Plots are generated for each sample using the CapCruncher package.

## Pipeline Output

The `capcruncher_output/results` directory contains the following files:

### General

- `design_matrix.tsv`: the design matrix used or generated by the pipeline. This
  file is used for the differential analysis and plot generation.

### Statistics and QC

- `capcruncher_report.html`: the main report of the pipeline. It contains the
  results of the analysis and the figures.

- `full_qc_report.html`: the full QC report of the pipeline. It contains the
  results of the QC analysis and the figures.

### Individual Samples

- `<SAMPLE_NAME>`: the results of the analysis for each sample. The directory contains the following files:

  - `bigwigs`: the bigwig files of the sample. The files are stored in the
      `raw` and `norm` directories. The `raw` directory contains the raw
      bigwig files. The `norm` directory contains the normalized bigwig files.
      **Note**: Only for Capture-C and Tri-C

  - `<SAMPLE_NAME>.bam`: the alignment file of the sample.

  - `<SAMPLE_NAME>.bam.bai`: the index of the alignment file of the sample.

  - `<SAMPLE_NAME>.hdf5`: Cooler formated groups within the HDF5 file per viewpoint. See [Cooler documentation](https://cooler.readthedocs.io/en/latest/) for more information on the format.

  - `<SAMPLE_NAME>.parquet`: This contains all of the data from the sample in a tabular format that can be accessed by any software that reads parquet files e.g. `pandas`, `polars`, [Arrow for R](https://github.com/apache/arrow/tree/main/r).

### Comparisons and Differential Analysis

!!! note
    The comparisons and differential analysis are only performed if two or more replicates are provided.
    Currently only Capture-C and Tri-C data are supported.

- `comparisons`: the results of the comparisons between the different
  conditions. The results are stored in the `bigwigs`
  directories

- `differential`: the results of the differential analysis

### Visualisation

- `figures`: Plots generated by the pipeline. for each viewpoint at the coordinates provided in the configuration file.
    This also generates PlotNado TOML templates that can be used with the `capcruncher plot` command. For more advanced or customisable plots, edit these templates directly or build figures with [PlotNado](https://alsmith151.github.io/plotnado/), which provides the underlying track options, template format, and Python API used by CapCruncher.

- UCSC Hub (if enabled in the configuration file): the UCSC Genome Browser track hub
  of the pipeline. It contains the bigwig files of the samples and the
  comparisons between samples (**Note**: Only for Capture-C and Tri-C)

!!! tip
        Unlike regular cooler.hdf5 files, there are multiple groups per file, one per viewpoint. You can use any cooler supported tool but you will need to specify the group name. For example, to get the matrix for the viewpoint `Slc25A37` you can use the following command:

        ``` bash
        cooler dump -t pixels -o Slc25A37.tsv <SAMPLE_NAME>.hdf5::Slc25A37`

        ```
