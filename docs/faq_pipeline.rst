Pipeline FAQ
############

Interruptions to the pipeline
=============================

If for any reason the pipeline interrupted, simply delete any malformed file or
directory and re-run the pipeline. The pipeline will then continue with the analysis
without needing to start from the beginning.

Leftover ctmpXXX.sh files can be safely deleted if the pipeline has stopped running::

    rm -f *.sh*

To restart the pipeline from the beginning for any reason, simply delete all
processed files and re-run the pipeline::

    rm -rf capcruncher_preprocessing/ capcruncher_analysis/ capcruncher_compare/ statistics/ pipeline.log
    capcruncher pipeline make


DRMAA error
===========

A common error when running the pipeline is:

::

    GLOBAL_SESSION = drmaa.Session()
    NameError: name 'drmaa' is not defined

This error occurrs because you are not connected to the cluster.
See :ref:`Pre-Installation recommendations`  for how to add DRMAA_LIBRARY_PATH to your bashrc.
Alternatively, if you do not want to use a cluster you can run the pipeline in :ref:`local mode <Step 4 - Running the pipeline>`. 

KeyError: 'missing parameter accessed'
======================================

This occurs if one of the keys in config.yml has been altered or deleted.
Either find the source of the error or get a new version of the yaml file and fill it out.

ValueError: not enough values to unpack (expected 2, got 1)
===========================================================

The fastq files are not paired correctly. Please ensure that fastq files are in
the format::

    NAME-OF-SAMPLE_(R)1.fastq.gz
    NAME-OF-SAMPLE_(R)2.fastq.gz
