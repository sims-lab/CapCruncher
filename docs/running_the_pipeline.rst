
Pipeline
########

The main feature of ccanalyser is the end-to-end data processing pipeline. 
The pipeline has been written using the `cgat-core pipelining library <https://github.com/cgat-developers/cgat-core>`_ 
and the following diagram illustrates the steps performed by the pipeline:

.. image:: images/pipeline_flow.svg
    :width: 100%
    :alt: Pipeline flow diagram


This section provides further details on how to run the pipeline. In essence
the pipeline requires a working directory with correctly named fastq files
and a :ref:`config.yml <Step 2 - Edit a copy of config.yml>` file that provides
the pipeline configuration.  



Step 1 - Create a working directory
===================================

To run the pipeline you will need to create a :term:`working directory`
for the pipeline run:

::

   mkdir RS411_EPZ5676/
   cd RS411_EPZ5676/

The pipeline will be executed here and all files will be generated
in this directory.

Step 2 - Edit a copy of config.yml
==================================

The configuration file `config.yml <https://github.com/sims-lab/capture-c/blob/master/config.yml>`_ enables 
parameterisation of the pipeline run with user specific settings. Furthermore,
it also provides paths to essential files for the pipeline run e.g. bowtie2 indicies.
The paths supplied do not have to be in the same directory as the pipeline but it is
recomended to copy the capture oligos used to the :term:`working directory`.

.. warning::

    The yaml file must be named **config.yml** for the pipeline to recognise it and run correctly.

.. code :: yaml

    ###################################
    # Essential configuration options #
    ###################################

    analysis:
        
        # Method to use for data analysis, choose from capture | tri | tiled
        method: capture 
        
        # Path to capture oligos used for analysis. This must be in bed format and named. 
        capture_oligos: capture-c_oligos.bed
        
        # Restriction enzyme name or recognition site
        restriction_enzyme: dpnii
        
        # Number of basepairs to exclude around each capture probe to remove re-ligations
        reporter_exclusion_zone: 1000

        # Genomic window size(s) to use for binning restriction fragment interaction counts
        # into even genomic windows.
        # Only bin sizes that are present here will be allowed to be used for heatmap generation.
        bin_size: 2500, 5000

    genome:
        
        # Name of genome. UCSC genome names are prefered.
        name: mm9
        
        # Path to fasta file containing entire genome sequence separated by chromosome. 
        fasta: /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa

        # Path to indicies for the specified aligner (default = bowtie2)
        aligner_index: /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome

        # Chromosome sizes for genome. Will be determined automatically from the genome name if this is a UCSC genome.
        chrom_sizes: 

    #############################
    # Essential cluster options #
    #############################

    pipeline:

        # Cluster manager used. e.g. SLURM, SunGrid Engine, etc.
        cluster_queue_manager: slurm

        # Name of queue to use. This will be cluster specific.
        cluster_queue: batch
        
        # Maximum number of cores to use per task. (~ 4 works best)
        n_cores: 4

        # Maximum memory avalible per job.
        # Some tasks are quite demanding on memory for a deeply sequenced experiment (32 G recomended)
        memory: 32G


    ###################################
    # Optional configuration options #
    ###################################

    align:
        
        # Aligner to use. Both bowtie and bowtie2 are supported but bowtie2 is prefered.
        aligner: Bowtie2

        # Flag to specify index for the aligner. Leave blank if this is not present i.e. bowtie
        index_flag: -x 

        # Aligner specific options. (Take care as this will ignore pipeline_n_cores) 
        options: -p 6 --very-sensitive

    analysis_optional:
        
        # Path to blacklisted regions bed file. Must be named. Can supply any regions to be removed.
        blacklist:

    deduplication:

        # Turns on initial removal of identical reads
        pre-dedup: True 

    hub:
        
        # Determines if hub is created or not. 
        create: False
        
        # Url/IP of server to host bigWigs
        url: http://userweb.molbiol.ox.ac.uk/
        
        # Location of publically accessible location on the server
        dir: /public/asmith/capture-c/test 
        
        # Name for the hub (UCSC required)
        name: capturec_test 
        
        # Short hub name (UCSC required)
        short: capturec new pipeline 
        
        # Long hub name (UCSC required)
        long: capturec processed with the new python pipeline 
        
        # Email address (UCSC required)
        email: alastair.smith@ndcls.ox.ac.uk 
        
        # Colours to use for bigWig tracks. Leave blank for random. Colours are cycled if there are more tracks than colours.
        colors: red blue green yellow purple 

    normalisation:
        
        # Scaling factor for normalisation. (not currently used but will be added shortly)
        scale_factor: 1000000 

    plot:

        # Path to a bed file containing coordinates for plotting a heatmap of reporters.
        # Must be named and the interval name must contain the probe name to be plotted.
        coordinates: plot_coords.bed

        # Plot output format (not currently used)
        format: png

        # Minimum value for plotting
        vmin: 0

        # Maximum value for plotting. Higher values are truncated to this value
        vmax: 1

        # Bin size(s) to use for plotting
        bin_size: 2500

        # Normalisation method to use for plot. Leave blank for default based on analysis method
        normalisation:

    trim:

        # Options passed to trim_galore
        options: --length 21 

    split:

        # Fastq files are split for parallel processing. Defines number of reads per fastq file (lower = more files to process)
        # For Tiled-C or highly enriched assays it is advised to reduce this to ~1e6 readpairs to speed up the analysis.
        # The pipeline will cope with a high number of reads but it will run slower. 
        n_reads: 2500000 


This yaml file can be edited using standard text editors e.g.:

::

    # To use gedit
    gedit config.yml

    # To use nano
    nano config.yml



Step 3 -  Copy or link fastq files into the :term:`working directory`
=====================================================================

The pipeline requires that fastq files are paired and in any of these formats:

.. note::
    
    Gziped files are handled appropriately without the need for extraction if .gz is
    present at the end of the file name. 

.. note::

    Multi-lane fastq files should be
    concatenated prior to running the pipeline otherwise multiple separate analyses will
    be performed.

* samplename_R1.fastq.gz
* samplename_1.fastq.gz
* samplename_R1.fastq
* samplename_1.fastq

All fastq files present in the directory will be processed by the pipeline in parallel and
original fastq files will not be modified. If new fastq files are added to a pre-run pipeline,
only the new files will be processed.



Copy:

::

    cp PATH_TO_FASTQ/example_R1.fastq.gz.

Symlink example:

Be sure to use the absolute path for symlinks

::

    ln -s ABSOLUTE_PATH_TO_FASTQ/example_R1.fastq.gz


Step 4 - Running the pipeline
=============================

After copying/linking fastq files into the working directory and configuring
config.yml for the current experiment, the pipeline can be run with:

::

    ccanalyser pipeline


There are several options to visualise which tasks will be performed by the pipeline
before running. 

The tasks to be performed can be examined with:

::
    
    # Shows the tasks to be performed
    ccanalyser pipeline show 

    # Plots a directed graph using graphviz
    ccanalyser pipeline plot

If you are happy with the tasks to be performed, the full pipeline run can be started with:

::

    # If using all default settings and using a cluster
    ccanalyser pipeline make

    # If not using a cluster, run in local mode.
    ccanalyser pipeline make --local -p 4

    # Avoiding disconnects
    nohup ccanalyser pipeline make &


See `cgatcore <https://cgat-core.readthedocs.io/en/latest/getting_started/Examples.html>`_ for additional
information.



Step 5 - Running the pipeline to a specified stage
==================================================

There are currently multiple stopping points built into the pipeline at key stages. These are:

* fastq_preprocessing - Stops after in silico digestion of fastq files.
* pre_annotation - Stops before aligned slices are ready to be annotated.
* post_annotation - Stops after aligned slices have been annotated.
* post_ccanalyser_analysis - Stops after reporters have been identified and duplicate filtered.
* full - Run the pipline until all required tasks are complete

To run the pipeline until one of these stopping points, use:

::

    # Run until TASK_NAME step
    ccanalyser pipeline make TASK_NAME



