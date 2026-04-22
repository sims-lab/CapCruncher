# Installation

!!! warning
    CapCruncher is currently only availible for linux. MacOS support is planned for the future.

## Setup

CapCruncher requires Python 3.12 or newer.

For end users, it is still highly recommended to install CapCruncher in a conda-compatible environment because the pipeline depends on several native bioinformatics tools. If you do not have conda installed, see the detailed [conda installation section](#detailed-conda-installation-instructions).

## Dependencies

There are two main ways to obtain the dependencies required to run CapCruncher:

### Install all dependencies using conda

#### Direct Installation

The easiest way to install these dependencies is to use conda. Run the following command to install CapCruncher and all dependencies:

```{bash}
mamba install -c bioconda capcruncher
```

!!! warning
    The latest version of CapCruncher is not yet available on conda. Please install the latest version from PyPI using the command below.


#### Two-step installation using conda and pip

Alternatively, create a new conda environment from the checked-in environment definition and install CapCruncher with pip (currently the recommended method):


```{bash}
wget https://raw.githubusercontent.com/sims-lab/CapCruncher/master/environment.yml
conda env create -f environment.yml
conda activate cc
python -m pip install --upgrade pip

# Install CapCruncher using pip
python -m pip install capcruncher

# Optional - highly recommended to install the optional dependencies
# Installs the plotting and analysis extras defined by the package.
python -m pip install 'capcruncher[full]'
```

The bundled [environment.yml](https://raw.githubusercontent.com/sims-lab/CapCruncher/master/environment.yml) now creates a Python 3.12 environment and installs the non-Python tools needed by the pipeline.

### Editable development install with uv

For local development, `uv` is the shortest way to get an isolated editable install on Python 3.12+:

```{bash}
git clone https://github.com/sims-lab/CapCruncher.git
cd CapCruncher

uv python install 3.12
uv venv --python 3.12 .venv
source .venv/bin/activate
uv pip install -e '.[full]'

# Cheap smoke check for the editable install
capcruncher --help
```

Use this path when you are changing Python code or running the test suite. Use the conda environment above when you need the full set of external command-line tools that the pipeline shells out to.

### Install CapCruncher in a minimal conda environment and use singularity to run the pipeline

!!! note
    Singularity must be installed on your system to use this method. See the [pipeline guide](pipeline.md) for more information. The pipeline will only function correctly if using the --use-singularity option. This is because the pipeline uses singularity containers to run the pipeline steps. See the [pipeline guide](pipeline.md) for more information.


Create a minimal conda environment and install CapCruncher using pip:

```{bash}
mamba create -n cc python=3.12
conda activate cc
python -m pip install --upgrade pip

# Install the package and optional extras from PyPI
python -m pip install 'capcruncher[full]'
```

If you are working from a local checkout instead of PyPI, replace the last command with:

```{bash}
# Optional - highly recommended to install the optional dependencies
python -m pip install -e '.[full]'
```


## Manual Installation (Not Recommended)

### Install Dependencies

See the dependencies in the [environment.yml](https://raw.githubusercontent.com/sims-lab/CapCruncher/master/environment.yml) and [requirements.txt](https://raw.githubusercontent.com/sims-lab/CapCruncher/master/requirements.txt) files. All dependencies can be installed using conda or pip.

### Install CapCruncher from GitHub

Clone the repository and install CapCruncher using pip:

```{bash}
git clone https://github.com/sims-lab/CapCruncher.git
cd CapCruncher
python -m pip install --upgrade pip
python -m pip install .

# Optional - highly recommended to install the optional dependencies
python -m pip install '.[full]'
```


## Detailed Conda Installation Instructions

Download and install MambaForge by following the [Miniforge installation instructions](https://github.com/conda-forge/miniforge#mambaforge) for your system. You will typically need the x86_64 installer on most Linux systems.

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
