# Installation

!!! warning
    CapCruncher is currently only available for Linux. MacOS support is planned for the future.

## Setup

It is highly recommended to install CapCruncher in a conda environment. If you do not have conda installed, see the detailed [conda installation section](#detailed-conda-installation).

## Dependencies

There are two main ways to obtain the dependencies required to run CapCruncher:

### Use the Docker image

CapCruncher publishes a Docker image containing the CLI, pipeline runtime, native tools, and reporting dependencies:

```{bash}
docker pull ghcr.io/sims-lab/capcruncher:latest
docker run --rm ghcr.io/sims-lab/capcruncher:latest --help
```

To run a pipeline from a directory containing `capcruncher_config.yml` and input FASTQs:

```{bash}
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -v "$PWD":/work \
  -w /work \
  ghcr.io/sims-lab/capcruncher:latest \
  pipeline --cores 8
```

See the [Docker guide](docker.md) for path mounting, local image builds, and examples for other CLI commands.

Docker is intended for local workstation and CI usage. On HPC systems, use
Apptainer; Docker daemons are generally not available or permitted on shared
clusters.

### Install all dependencies using conda

#### Direct Installation

The easiest way to install these dependencies is to use conda. Run the following command to install CapCruncher and all dependencies:

```{bash}
mamba install -c bioconda capcruncher
```

!!! warning
    The latest version of CapCruncher is not yet available on conda. Please install the latest version from PyPI using the command below.


#### Two-step installation using conda and pip

Alternatively, create a new conda environment and install CapCruncher using pip (currently the recommended method):


```{bash}
wget https://raw.githubusercontent.com/sims-lab/CapCruncher/master/environment.yml
conda env create -f environment.yml
conda activate cc

# Install CapCruncher using pip
pip install capcruncher

# Install the full dependency set used by the pipeline and CLI.
pip install capcruncher[full]
```

### Install CapCruncher in a minimal conda environment and use Apptainer to run the pipeline

!!! note
    Apptainer is the supported container runtime for HPC usage. See the [pipeline guide](pipeline.md) for more information. Use the bundled `capcruncher-local-apptainer` or `capcruncher-slurm-apptainer` presets to run workflow steps in containers.


Create a minimal conda environment and install CapCruncher using pip:

```{bash}
mamba create -n cc "python>=3.12"
conda activate cc
pip install capcruncher[full]
```

Install the editable Snakemake profiles after installation:

```{bash}
capcruncher pipeline-init
```

This writes profiles to `${XDG_CONFIG_HOME:-~/.config}/snakemake`, Snakemake's
standard user profile directory on Linux. See the [cluster setup guide](cluster_config.md)
for editing profiles and refreshing them with `capcruncher pipeline-init --force`.

You can also run the CapCruncher container directly with Apptainer:

```{bash}
apptainer exec docker://ghcr.io/sims-lab/capcruncher:latest capcruncher --help
```


## Manual Installation (Not Recommended)

### Install Dependencies

See the dependencies in the [environment.yml](https://raw.githubusercontent.com/sims-lab/CapCruncher/master/environment.yml) and [requirements.txt](https://raw.githubusercontent.com/sims-lab/CapCruncher/master/requirements.txt) files. All dependencies can be installed using conda or pip.

### Install CapCruncher from GitHub

Clone the repository and install CapCruncher using pip:

```{bash}
git clone https://github.com/sims-lab/CapCruncher.git
cd CapCruncher
pip install .

# Install the full dependency set used by the pipeline and CLI.
pip install .[full]
```


## Detailed Conda Installation Instructions

Download and install MambaForge from [here](https://github.com/conda-forge/miniforge#mambaforge) for your system (You will typically need the x86_64 version for most Linux systems).

### Download and run the installer for your system (only Linux is supported at the moment)

```{bash}
# Download the installer for your system
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh

# Allow the installer to be executed
chmod +x Mambaforge-Linux-x86_64.sh

# Run the installer
./Mambaforge-Linux-x86_64.sh
```

Follow the instructions to install MambaForge. It is advised to install MambaForge in a location with a reasonable amount of free space (>2GB) as it will be used to install all dependencies for CapCruncher.

### Initialise MambaForge in your shell

```{bash}
conda init bash
```

### Refresh your shell

```{bash}
source ~/.bashrc
```


### Setup conda channels

```{bash}
conda config --set channel_priority strict
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Now the installation installation of CapCruncher can be completed using the instructions [above](#dependencies).
