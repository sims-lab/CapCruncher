# Installation

!!! warning
    CapCruncher is currently only availible for linux. MacOS support is planned for the future.

## Setup

It is highly recommended to install CapCruncher in a conda environment. If you do not have conda installed, see the detailed [conda installation section](#detailed-conda-installation).

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

Alternatively, create a new conda environment and install CapCruncher using pip (currently the recommended method):


```{bash}
wget https://raw.githubusercontent.com/sims-lab/CapCruncher/master/environment.yml
conda env create -f environment.yml
conda activate cc

# Install CapCruncher using pip
pip install capcruncher

# Optional - highly recommended to install the optional dependencies
pip install capcruncher[stats,plotting,experimental]
```

### Install CapCruncher in a minimal conda environment and use singularity to run the pipeline

!!! note
    Singularity must be installed on your system to use this method. See the [pipeline guide](pipeline.md) for more information. The pipeline will only function correctly if using the --use-singularity option. This is because the pipeline uses singularity containers to run the pipeline steps. See the [pipeline guide](pipeline.md) for more information.


Create a minimal conda environment and install CapCruncher using pip:

```{bash}
mamba create -n cc "python>=3.10"
conda activate cc
# Optional - highly recommended to install the optional dependencies
pip install capcruncher[stats,plotting,experimental]
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

# Optional - highly recommended to install the optional dependencies
pip install .[stats,plotting,experimental]
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
