<p align="center">
  <img src="docs/img/capcruncher_logo.png" alt="CapCruncher" width="760">
</p>

<p align="center">
  <a href="https://sims-lab.github.io/CapCruncher/"><img alt="Documentation" src="https://github.com/sims-lab/CapCruncher/actions/workflows/docs.yml/badge.svg"></a>
  <a href="https://codecov.io/gh/sims-lab/CapCruncher"><img alt="Codecov" src="https://codecov.io/gh/sims-lab/CapCruncher/graph/badge.svg?token=RHIGNMGX09"></a>
  <a href="https://anaconda.org/bioconda/capcruncher"><img alt="Bioconda" src="https://anaconda.org/bioconda/capcruncher/badges/version.svg"></a>
  <a href="https://anaconda.org/bioconda/capcruncher"><img alt="License" src="https://anaconda.org/bioconda/capcruncher/badges/license.svg"></a>
  <a href="https://zenodo.org/badge/latestdoi/224631087"><img alt="DOI" src="https://zenodo.org/badge/224631087.svg"></a>
  <a href="https://pepy.tech/project/capcruncher"><img alt="Downloads" src="https://pepy.tech/badge/capcruncher"></a>
</p>

CapCruncher is a command-line toolkit and Snakemake workflow for processing
Capture-C, Tri-C and Tiled-C sequencing data. It provides capture-aware
filtering, contact matrix generation, plotting, and UCSC track hub output for
workstation, container, and HPC runs.

See the [documentation](https://sims-lab.github.io/CapCruncher/) for the full
installation, configuration, and pipeline guides.

## What It Does

- Processes raw Capture-C, Tri-C and Tiled-C FASTQs into filtered contacts.
- Builds contact matrices and viewpoint-centric outputs for downstream analysis.
- Generates PlotNado-compatible plot templates and rendered figures.
- Produces UCSC genome browser track hubs.
- Runs locally, in Docker, or on HPC systems through Snakemake 9 presets.

## Installation

The recommended native install is a Conda/Mamba environment from Bioconda:

```bash
mamba create -n capcruncher -c conda-forge -c bioconda capcruncher
conda activate capcruncher
capcruncher --help
```

On HPC systems, use Apptainer — it runs rootless and integrates with Slurm.
Apptainer pulls the image automatically from the registry:

```bash
apptainer exec docker://ghcr.io/sims-lab/capcruncher:latest capcruncher --help
```

On workstations, the Docker image bundles the native bioinformatics tools with
the Python runtime:

```bash
docker run --rm ghcr.io/sims-lab/capcruncher:latest --help
```

CapCruncher is also published to PyPI for Python-side CLI/API usage, or for
systems where the native tools are already installed:

```bash
pip install capcruncher
```

Pixi is used for development and CI reproducibility, but it is not the standard
end-user install route. CapCruncher targets Linux execution. macOS users can run
the Linux container via Docker Desktop or Colima. For fallback native installs,
Apptainer profiles, and development setup, see the
[installation guide](docs/installation.md).

## Quick Start

Create a pipeline configuration:

```bash
capcruncher pipeline config
```

Install the bundled Snakemake profiles:

```bash
capcruncher pipeline init
```

Run from the directory containing `capcruncher_config.yml` and your FASTQ files:

```bash
capcruncher pipeline run --cores 8 --preset local
```

For an HPC run with Apptainer-backed jobs:

```bash
capcruncher pipeline run --jobs 50 --preset slurm-apptainer
```

For Docker-based workstation usage:

```bash
docker run --rm -it -v "$PWD":/work -w /work \
  ghcr.io/sims-lab/capcruncher:latest \
  pipeline --cores 8
```

Long cluster runs should be launched inside `tmux`, `screen`, or your scheduler's
normal job submission wrapper so they survive terminal disconnects.

## CLI

List available commands:

```bash
capcruncher --help
```

Inspect a command:

```bash
capcruncher <command> --help
```

See the [CLI reference](https://sims-lab.github.io/CapCruncher/cli/) and
[pipeline guide](https://sims-lab.github.io/CapCruncher/pipeline/) for detailed
usage.

## Development

This repository currently targets Python 3.12+ and Snakemake 9. The local
validation environment used by maintainers is the `cc` conda environment:

```bash
conda run -n cc pytest tests/test_cli.py -q
conda run -n cc pytest tests/test_workflow_scripts.py tests/test_plotting.py -q
```

Please report bugs and feature requests through the
[issue tracker](https://github.com/sims-lab/CapCruncher/issues/new/choose).
