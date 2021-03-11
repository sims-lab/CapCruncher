====================
Running the pipeline
====================

This section provides details on how to run the pipeline.

Pipeline set-up
======================

**Step 1**: To run the pipeline you will need to create a :term:`working directory`
and enter it. For example:

.. code-block::bash
   mkdir cc_pipeline/
   cd cc_pipeline/

The pipeline will be executed here and all files will be generated
in this directory.

**Step 2**: Edit a copy of config.yml:

The YAML file allows the pipeline to be parameterised for different situations.
Furthermore the pipeline requires a number of additional files in order to run.
e.g. Paths to whole genome fasta files, bowtie2 indicies, capture_oligos.

This can then be edited using standard text editors e.g.:

.. code-block::bash
    # To use gedit
    gedit config.yml

    # To use nano
    nano config.yml

**Note** - currently the yaml file must be named config.yml for the
pipeline to recognise it and run correctly.

**Step 3**: Copy or link fastq files into the :term:`working directory`:

The pipeline requires that fastq files are paired and in any of these formats:

    * samplename_R1.fastq.gz
    * samplename_1.fastq.gz
    * samplename_R1.fastq
    * samplename_1.fastq

All fastq files present in the directory will be processed by the pipeline.
Original fastq files will not be modified and will not be re-run if processing
has already been completed. If new fastq files are added to a pre-run pipeline,
only the new files will be processed.

Gziped files are handled appropriately without the need for extraction if .gz is
present at the end of the file name. Multi-lane fastq files should be
concatenated prior to running the pipeline otherwise multiple separate analyses will
be performed.

Copy:

.. code-block::bash
    cp PATH_TO_FASTQ/example_R1.fastq.gz.

Symlink example:
**Note** - Be sure to use the absolute path for symlinks

.. code-block::bash
    ln -s ABSOLUTE_PATH_TO_FASTQ/example_R1.fastq.gz


Running the pipeline
=====================

When all components are present, the tasks to be run can be:

Shown with:

.. code-block::bash

    ccanalyser pipeline show

Plotted in inkscape with:

.. code-block::bash

    ccanalyser pipeline plot

Run with:

.. code-block::bash

    # If using all default settings and using a cluster
    ccanalyser pipeline make

    # If not using a cluster
    ccanalyser pipeline make --local -p 4

Troubleshooting
===============

**Interruptions to the pipeline**

If for any reason the pipeline interrupted, simply delete any malformed file or
directory and re-run the pipeline. The pipeline will then continue with the analysis
without needing to start from the beginning.

Leftover ctmpXXX.sh files can be safely deleted if the pipeline has stopped running.

.. code-block::bash
    rm -f *.sh*

To restart the pipeline from the beginning for whatever reason, simply delete all
processed files and re-run the pipeline:

.. code-block::bash
    rm -rf pre-ccanalysis/ ccanalysis/ capture_compare/ statistics/ pipeline.log
    ccanalyser pipeline make


**DRMAA error**

One of the most common errors when running the pipeline is:

.. code-block::bash

    GLOBAL_SESSION = drmaa.Session()
    NameError: name 'drmaa' is not defined

This error occurrs because you are not connected to the cluster.
See :ref:`Pre-Installation recommendations`  for how to add DRMAA_LIBRARY_PATH to your bashrc.
Alternatively, if you do not want to use a cluster (or can't) you can run the
pipeline in local mode see above.

**KeyError: 'missing parameter accessed'**

This occurs if one of the keys in capturec_pipeline.yml has been altered or deleted.
Either find the source of the error or get a new version of the yaml file and fill it out.

**ValueError: not enough values to unpack (expected 2, got 1)**

The fastq files are not paired correctly. Please ensure that fastq files are in
the format:

    | NAME-OF-SAMPLE_(R)1.fastq.gz
    | NAME-OF-SAMPLE_(R)2.fastq.gz
