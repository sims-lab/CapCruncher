===================================
Pipeline
===================================

This section provides details on how to run the pipeline.

Pipeline set-up
======================

**Step 1**: To run the pipeline you will need to create a :term:`working directory`
and enter it. For example:

::

   mkdir cc_pipeline/
   cd cc_pipeline/

The pipeline will be executed here and all files will be generated
in this directory.

**Step 2**: Edit a copy of config.yml:

The YAML file allows the pipeline to be parameterised for different situations.
Furthermore the pipeline requires a number of additional files in order to run.
e.g. Paths to whole genome fasta files, bowtie2 indicies, capture_oligos.

This can then be edited using standard text editors e.g.:

::

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

::

    cp PATH_TO_FASTQ/example_R1.fastq.gz.

Symlink example:
**Note** - Be sure to use the absolute path for symlinks

::

    ln -s ABSOLUTE_PATH_TO_FASTQ/example_R1.fastq.gz


Running the pipeline
=====================

When all components are present, the tasks to be run can be:

Shown with:

::

    ccanalyser pipeline show

Plotted in inkscape with:

::

    ccanalyser pipeline plot

Run with:

::

    # If using all default settings and using a cluster
    ccanalyser pipeline make

    # If not using a cluster
    ccanalyser pipeline make --local -p 4

