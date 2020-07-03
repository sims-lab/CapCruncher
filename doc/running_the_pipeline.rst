====================
Running the pipeline
====================

This section provides details on how to run the pipeline.

Pipeline set-up
======================

**Step 1**: To run the pipeline you will need to create a :term:`working directory`
and enter it. For example::

   mkdir cc_pipeline/
   cd cc_pipeline/

The pipeline will be executed here and all files will be generated
in this directory.

The original fastq files will not be modified.

**Step 2**: Edit a copy of capturec_pipeline.yml:

The YAML file allows the pipeline to be parameterised for different situations.
Furthermore the pipeline requires a number of additional files in order to run.
See capturec_pipeline.yml for further details.::

   cp PATH_TO_REPO/capture-c/capturec_pipeline.yml .
   gedit capturec_pipeline.yml

**Step 3**: Copy or link fastq files into the :term:`working directory`:

Copy example
::
    cp PATH_TO_FASTQ/example_R1.fastq.gz.

Symlink example
::
    ln -s ABSOLUTE_PATH_TO_FASTQ/example_R1.fastq.gz

The pipeline requires that fastq files are paired and in any of these formats:

    * samplename_R1.fastq.gz
    * samplename_1.fastq.gz
    * samplename_R1.fastq
    * samplename_1.fastq

Gziped files are handled appropriately without the need for extraction if .gz is
present at the end of the file name. Multi-lane fastq files should be
concatenated prior to running the pipeline.

All fastq files present in the directory will be processed by the pipeline. They
will not be re-run if processing has already been completed i.e. if new fastq
files are added to a pre-run pipeline, only the new files will be processed.
