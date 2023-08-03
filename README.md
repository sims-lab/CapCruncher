# CapCruncher

[![Documentation Status](https://readthedocs.org/projects/capcruncher/badge/?version=latest)](https://sims-lab.github.io/CapCruncher/)
[![codecov](https://codecov.io/gh/sims-lab/CapCruncher/branch/master/graph/badge.svg?token=RHIGNMGX09)](https://codecov.io/gh/sims-lab/CapCruncher)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/capcruncher/badges/version.svg)](https://anaconda.org/bioconda/capcruncher)
[![Anaconda-Server Badge License](https://anaconda.org/bioconda/capcruncher/badges/license.svg)](https://anaconda.org/bioconda/capcruncher)
[![DOI](https://zenodo.org/badge/224631087.svg)](https://zenodo.org/badge/latestdoi/224631087)
[![Downloads](https://pepy.tech/badge/capcruncher)](https://pepy.tech/project/capcruncher)

![CapCruncher Logo](docs/img/capcruncher_logo.png)

CapCruncher is a package explicitly designed for processing Capture-C, Tri-C and Tiled-C data. Unlike other pipelines that are designed to process Hi-C or Capture-HiC data, the filtering steps in CapCruncher are specifically optimised for these datasets.

The package consists of a configurable data processing pipeline and a supporting command line interface to enable fine-grained control over the analysis.

The pipeline is fast, robust and scales from a single workstation to a large HPC cluster. The pipeline is designed to be run on a HPC cluster and can be configured to use a variety of package management systems e.g. conda and singularity.

For more information, see the [documentation](https://sims-lab.github.io/CapCruncher/).

**Note:**
The current version of CapCruncher is in beta. Please report any issues you encounter to the [issue tracker](https://github.com/sims-lab/CapCruncher/issues/new/choose)


## Quick Start

### Installation

**Warning:**
CapCruncher is currently only availible for linux. MacOS support is planned for the future.

CapCruncher is available on conda and PyPI. To install the latest version, run:

It is highly recommended to install CapCruncher in a conda environment. If you do not have conda installed, please follow the instructions [here](https://github.com/conda-forge/miniforge#mambaforge) to install mambaforge.

``` bash
pip install capcruncher
```

or

``` bash
mamba install -c bioconda capcruncher
```

See the [installation guide](installation.md) for more detailed instructions.

### Usage

CapCruncher commands are run using the `capcruncher` command. To see a list of available commands, run:

``` bash
capcruncher --help
```

To see a list of available options for a command, run:

``` bash
capcruncher <command> --help
```

See the [usage guide](usage.md) for more detailed instructions.

### Pipeline

The CapCruncher pipeline handles the processing of raw data from the sequencer to the generation of a contact matrix, generation of plots and production of a UCSC genome browser track hub.

See the [pipeline guide](pipeline.md) for more detailed instructions including how to configure the pipeline to run on HPC clusters and using various package management systems e.g. conda and singularity.

#### Pipeline Configuration

The pipeline is configured using a YAML file. It is strongly recommended to use the `capcruncher pipeline-config` command to generate a template configuration file. This command will generate a template configuration file with all available options and descriptions of each option.

``` bash
capcruncher pipeline-config --help
```

#### Running the pipeline

The pipeline is run using the `capcruncher pipeline` command. Ensure that you have a configuration file and the fastq files to process are in the current working directory.

``` bash
capcruncher pipeline --cores <NUMBER OF CORES TO USE>
```
