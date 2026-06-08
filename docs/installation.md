# Installation

!!! warning
    CapCruncher targets Linux execution. macOS users should run the Linux
    container through Docker Desktop, Colima, or Apptainer.

## Which install method should I use?

| Situation | Recommended method |
| --- | --- |
| HPC cluster | [Apptainer](#highly-recommended-containers) — pull once, no root required |
| HPC cluster / Linux workstation | [Bioconda / Mamba](#recommended-native-install) — single command |
| macOS workstation | [Docker](#highly-recommended-containers) — pipeline tools unavailable natively |
| Bioconda lags latest release | [conda + uv fallback](#fallback-native-install) |
| Python analysis only (no pipeline) | [pip](#python-only-install) |
| Development / contributing | [Pixi](#developer-install) |

## Recommended Native Install

The easiest native install is a Conda/Mamba environment from Bioconda:

```bash
mamba create -n capcruncher -c conda-forge -c bioconda capcruncher
conda activate capcruncher
capcruncher --help
```

This route installs CapCruncher with the native command-line tools required by
the pipeline, including aligners, FASTQ/QC tools, samtools, and UCSC utilities.
Use strict channel priority if you manage channels globally:

```bash
conda config --set channel_priority strict
```

If the Bioconda package is behind the latest PyPI release, use the fallback
environment below.

## Highly Recommended: Containers

Containers avoid native dependency conflicts and are highly recommended for
reproducible pipeline runs.

### Apptainer (HPC)

Apptainer runs rootless on most HPC clusters. For cluster-scale runs, install
the bundled Snakemake profiles and use the Apptainer-backed preset — Apptainer
will pull and cache the image automatically on the head node:

```bash
capcruncher pipeline init
capcruncher pipeline run --preset capcruncher-slurm-apptainer --jobs 50
```

On clusters where compute nodes lack internet access, pull the image to a `.sif`
file on the head node first, then point the config at the local file:

```bash
# Pull once on the head node (requires internet)
apptainer pull capcruncher.sif docker://ghcr.io/sims-lab/capcruncher:latest

# Set the local image path in capcruncher_config.yml
# execution:
#   container_image: /path/to/capcruncher.sif

capcruncher pipeline run --preset capcruncher-slurm-apptainer --jobs 50
```

For quick interactive use:

```bash
apptainer exec docker://ghcr.io/sims-lab/capcruncher:latest capcruncher --help
```

### Docker (workstations)

Use Docker on a local workstation:

```bash
docker pull ghcr.io/sims-lab/capcruncher:latest
docker run --rm ghcr.io/sims-lab/capcruncher:latest --help
```

Run a pipeline from a directory containing `capcruncher_config.yml` and input
FASTQs:

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -v "$PWD":/work \
  -w /work \
  ghcr.io/sims-lab/capcruncher:latest \
  pipeline run --cores 8
```

See the [Docker and Apptainer guide](docker.md) and
[cluster setup guide](cluster_config.md) for mount paths, profile editing, and
HPC examples.

## Fallback Native Install

If Bioconda is not current enough for your use case, create the maintained
environment and install the latest CapCruncher package from PyPI:

```bash
wget https://raw.githubusercontent.com/sims-lab/CapCruncher/master/environment.yml
mamba env create -f environment.yml
conda activate cc
uv pip install capcruncher
```

The environment file provides the native tools and Python runtime dependencies.
The final `uv pip install` step installs the current CapCruncher package.

## Python-Only Install

Use pip only when you already have the native tools installed, or when you only
need Python-side CLI/API commands:

```bash
pip install capcruncher
```

Optional Python features can be installed with extras:

```bash
pip install "capcruncher[plot,hub,differential,config,hpc]"
```

To install every Python-side optional feature:

```bash
pip install "capcruncher[all]"
```

Pure pip does not install native pipeline tools such as `bowtie2`, `samtools`,
`fastqc`, `flash2`, or UCSC command-line utilities.

## Developer Install

Pixi is the preferred developer and CI environment because it locks Python and
native dependencies reproducibly:

```bash
pixi install -e test
pixi run -e test pytest -q -m "not pipeline"
```

For editable development without Pixi, use the fallback native environment and
install the repository in editable mode:

```bash
mamba env create -f environment.yml
conda activate cc
pip install -e ".[plot,hub,differential,config,hpc]"
```

## Detailed Conda Setup

Install Miniforge or Mambaforge if you do not already have Conda/Mamba:

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
chmod +x Miniforge3-Linux-x86_64.sh
./Miniforge3-Linux-x86_64.sh
```

Initialise Conda for your shell and refresh it:

```bash
conda init bash
source ~/.bashrc
```

Then use the recommended native install command at the top of this page.
