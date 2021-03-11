====================
Installing ccanalyser
====================

Pre-Installation recommendations
===============================

1. Install conda if it has not been already using the `conda install instructions <https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent>`_.

2. If ccanalyser is **not** being installed through conda, first generate a new conda
environment using the yaml file in the `GitHub repo <https://github.com/sims-lab/capture-c/blob/master/capturec_conda_env.yml>`_:

.. code-block::bash
    conda env create -f ccanalyser_conda_env.yml
    conda activate capture-c

3. If you intend to use a cluster e.g. SunGrid engine/SLURM add the path to the DRMAA interface to your .bashrc:

.. code-block::bash
    # Access to the DRMAA library: https://en.wikipedia.org/wiki/DRMAA
    echo "export DRMAA_LIBRARY_PATH=/<full-path>/libdrmaa.so" >> ~/.bashrc

    # You can get this value from your configured environment:
    env | grep DRMAA_LIBRARY_PATH

    # or just look for the library:
    find / -name "*libdrmaa.so"


Installation
======================

The package can be installed in several ways:

.. 1. Install from conda:
.. .. code-block::bash
..     conda install ccanalyser

.. 2. Install from pypi:
.. .. code-block::bash
..     pip install ccanalyser

3. Install from GitHub:

.. code-block::bash
    git clone https://github.com/sims-lab/capture-c.git
    cd capture-c
    pip install .
