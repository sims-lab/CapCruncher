**********
Installing
**********

Recommendations prior to installation
#####################################

1) Install conda if it has not been already using the `conda install instructions <https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent>`_.

2) Ensure that conda channels are correctly set up:

::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

3) If capcruncher is **not** being installed through conda, first generate a new conda environment using the yaml file in the `GitHub repository <https://github.com/sims-lab/CapCruncher/blob/master/environment.yml>`_:

::

    wget https://raw.githubusercontent.com/sims-lab/CapCruncher/master/environment.yml
    conda env create -f environment.yml
    conda activate cc

4) (Highly Recommended) If you intend to use CapCruncher with a computational cluster e.g. SunGrid engine/SLURM add the path to the `DRMAA interface <https://en.wikipedia.org/wiki/DRMAA>`_ to your .bashrc:

::

    # You can get this value from your configured environment:
    env | grep DRMAA_LIBRARY_PATH

    # or just look for the library:
    find / -name "*libdrmaa.so"

    # To insert this path into your .bashrc
    echo "export DRMAA_LIBRARY_PATH=/<full-path>/libdrmaa.so" >> ~/.bashrc

    # An example for users of the CBRG computational cluster
    echo "export DRMAA_LIBRARY_PATH=/usr/lib64/libdrmaa.so" >> ~/.bashrc



Installation
############

The package can be installed in several ways:

1) Install from conda:

::

    conda install capcruncher

2) Install from pypi:

::

    pip install capcruncher

3) Install from GitHub:

::

    git clone https://github.com/sims-lab/CapCruncher.git
    cd CapCruncher
    pip install .


Installing optional packages
############################

To activate the plotting capabilities of the package additional dependencies are required.
In order to install these dependencies please run:

.. note::
    libcurl-dev must be present on your system in order to compile the plotting dependencies and
    may require administrator privileges to install.


::

    # Activate the CapCruncher conda environment
    conda activate cc

    # Install the plotting dependencies if using the pip/conda install method.
    pip install capcruncher[plotting]

    # Install the plotting dependencies if using the GitHub install method.
    cd PATH_TO_CAPCRUNCHER_REPOSITORY_CLONE/
    pip install .[plotting]
