# CapCruncher

[![Documentation](https://github.com/sims-lab/CapCruncher/actions/workflows/docs.yml/badge.svg?branch=master)](https://sims-lab.github.io/CapCruncher/)
[![Codecov](https://codecov.io/gh/sims-lab/CapCruncher/branch/master/graph/badge.svg?token=RHIGNMGX09)](https://codecov.io/gh/sims-lab/CapCruncher)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/capcruncher/badges/version.svg)](https://anaconda.org/bioconda/capcruncher)
[![Anaconda-Server Badge License](https://anaconda.org/bioconda/capcruncher/badges/license.svg)](https://anaconda.org/bioconda/capcruncher)
[![DOI](https://zenodo.org/badge/224631087.svg)](https://zenodo.org/badge/latestdoi/224631087)
[![Downloads](https://pepy.tech/badge/capcruncher)](https://pepy.tech/project/capcruncher)

![CapCruncher Logo](https://raw.githubusercontent.com/sims-lab/CapCruncher/master/docs/img/capcruncher_logo.png)

The CapCruncher package is designed to process Capture-C, Tri-C and Tiled-C data. Unlike other pipelines that are designed to process Hi-C or Capture-HiC data, the filtering steps in CapCruncher are specifically optimized for these datasets. The package consists of a configurable data processing pipeline and a supporting command line interface to enable fine-grained control over the analysis. The pipeline is fast, robust and scales from a single workstation to a large HPC cluster. It is designed to be run on an HPC cluster and can be configured to use a variety of package management systems, such as conda and singularity. For more information, see the [documentation](https://sims-lab.github.io/CapCruncher/).

> **Note:**
> The current version of CapCruncher is in beta. Please report any issues you encounter to the [issue tracker](https://github.com/sims-lab/CapCruncher/issues/new/choose)

## Quick Start

### Installation

> **Warning:**
>
> CapCruncher is currently only availible for linux with MacOS support planned in the future.

CapCruncher is available on conda and PyPI. To install the latest version, run:

``` bash
pip install capcruncher
```

or

``` bash
mamba install -c bioconda capcruncher
```

Please note that it is highly recommended to install CapCruncher in a conda environment. If you do not have conda installed, please follow the instructions [here](https://github.com/conda-forge/miniforge#mambaforge) to install mambaforge.

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

See the [CLI Reference](https://sims-lab.github.io/CapCruncher/cli/) for more detailed information regarding the various subcommands.

### Pipeline

The CapCruncher pipeline handles the processing of raw data from the sequencer to the generation of a contact matrix, generation of plots and production of a UCSC genome browser track hub. See the [pipeline guide](https://sims-lab.github.io/CapCruncher/pipeline/) for more detailed instructions including how to configure the pipeline to run on [HPC clusters](https://sims-lab.github.io/CapCruncher/pipeline/#hpc-cluster-usage-recommended-if-available) and use various package management systems [conda](https://sims-lab.github.io/CapCruncher/installation/#install-all-dependencies-using-conda) and [singularity](https://sims-lab.github.io/CapCruncher/pipeline/#singularity-usage-recommended-if-available).

#### Pipeline Configuration

The pipeline is configured using a YAML file but it is strongly recommended to use the `capcruncher pipeline-config` command to generate a tailored config file. To use the command, run:

``` bash
capcruncher pipeline-config
```

Simply follow the prompts to generate a config file. See the [pipeline configuration guide](https://sims-lab.github.io/CapCruncher/pipeline/#configuration-file) for more detailed instructions.

#### Running the pipeline

The pipeline is run using the `capcruncher pipeline` command. Ensure that you have a configuration file and the fastq files to process are in the current working directory.

``` bash
# Basic usage
capcruncher pipeline --cores <NUMBER OF CORES TO USE>

# Real example running the pipeline with 8 cores, using the slurm profile for running on a cluster with a SLURM workflow management system and using singularity for dependency management
capcruncher pipeline --cores 8 --profile slurm --use-singularity
```

> **Note:**
> In order to avoid disconnecting from the cluster, it is recommended to run the pipeline in a [tmux](https://linuxize.com/post/getting-started-with-tmux/)
> session. Alternatively, [nohup](https://linuxize.com/post/linux-nohup-command/) can be used to run the pipeline in the background. For example:
>
> ``` bash
> # tmux example
>tmux new -s capcruncher
> capcruncher pipeline --cores 8 --profile slurm --use-singularity
>
># nohup example
>nohup capcruncher pipeline --cores 8 --profile slurm --use-singularity &
>```

See the [pipeline guide](https://sims-lab.github.io/CapCruncher/pipeline/) for more detailed instructions.
