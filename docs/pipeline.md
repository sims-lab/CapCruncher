# Pipeline

The CapCruncher pipeline handles the processing of raw data from the sequencer to the generation of a contact matrix, generation of plots and production of a UCSC genome browser track hub.

This pipeline is based on the Snakemake workflow management system. Snakemake is a Python-based workflow management system that allows for the creation of reproducible and scalable data analyses. All elements of the workflow have been wrapped into the CapCruncher Python package. This allows for the pipeline to be run using the `capcruncher pipeline` command rather than having to run the pipeline using Snakemake directly.

## Pipeline Configuration

### Configuration File

The pipeline is configured using a YAML file. It is strongly recommended to use the `capcruncher pipeline-config` command to generate a template configuration file. This command will generate a template configuration file with all available options and descriptions of each option.

```bash
capcruncher pipeline-config
```

This utility will walk through the configuration options and generate a configuration file. It will generate a new directory <DATE>_<PROJECT NAME>_<ASSAY> and place the filled-out `capcruncher_config.yml` file in this directory.

For an example configuration file, see [here](examples/capcruncher_config.yml).

The configuration file can be edited manually if required e.g. to add a manually generated `design` file. Just ensure that the configuration file is valid YAML. A common error is to use tabs instead of spaces, this will cause the pipeline to fail while parsing the configuration file.

All options in the configuration file are documented within the file itself. Only the required options need to be filled out. The pipeline will use default values for all other options.

### Design File

The design file is a tab/comma/space-delimited file that contains the sample names and the metadata for each sample. This file is completely optional and only used for comparisons between Capture-C and Tri-C data. If it is not provided the pipeline will perform a basic sample name comparison to generate a basic design file. However, this will not be as accurate as a manually generated design file. The `design` file is a tab delimited file with the following columns:

- `sample`: The name of the FASTQ file (without the _R1.fastq.gz or _2.fastq.gz suffix)
- `condition`: The Group that the sample belongs to.

Provide the path to this file in the config file under the `design` key.


### Pipeline Steps


### All Assays

For all assays the pipeline consists of the following steps:

1. **Quality Control**: FastQC is used to perform quality control on the FASTQ files.
1. **Read Splitting**: The FASTQ files are split into parts of a user-defined size (default 1 million reads per part). This is done to allow for parallel processing of the FASTQ files and to reduce the memory requirements of the pipeline.
1. **Remove PCR Duplicates**: PCR duplicates are removed from the FASTQ files using the CapCruncher package to reduce the memory and CPU requirements of the pipeline.
1. **Read Trimming**: Trimming of the FASTQ files is performed using Trim Galore.
1. **Read Combining**: The trimmed FASTQ files are combined using FLASh to obtain the ligtion junctions.
1. **Read *in silico* Digestion**: The combined FASTQ files are digested *in silico* using the restriction enzyme or site specified in the configuration file.
1. **Read Alignment**: The digested FASTQ files are aligned to the reference genome using bowtie2.
1. **Alignment Annotation**: The aligned reads are annotated using the CapCruncher package.
1. **Alignment Filtering**: The aligned reads are filtered using the CapCruncher package.
1. **Alignment PCR Duplicate Removal**: PCR duplicates are removed from the aligned reads using the CapCruncher package.
1. **Contact Matrix Generation**: Contact matrices are generated using the CapCruncher package and stored in cooler (HDF5) format.
1. **Pipeline Statistics**: Statistics are generated for each sample using the CapCruncher package.
1. **Pipeline Plots**: Plots and `capcruncher plot` compatible templates (TOML format) are generated for each sample using the CapCruncher package.

### Capture-C and Tri-C

1. **Pileup Generation**: BigWig files are generated for each sample using the CapCruncher package.
1. **Pileup Normalisation**: The pileup files are normalised using the CapCruncher package.
1. **Pileup Comparison**: Interactions are compared between samples (if two or more replicates are provided) using the CapCruncher package.
1. **Differenital Interaction Analysis**: Differential interactions are identified between groups of samples using PyDESeq2.
1. **UCSC Genome Browser Track Hub**: A UCSC Genome Browser track hub is generated using the CapCruncher package.

### Tiled-C

1. **Contact Matrix Normalisation**: The contact matrices are normalised using the CapCruncher package with various third-party integrations.
1. **Plot Generation**: Plots are generated for each sample using the CapCruncher package.


## Running the pipeline

### Basic Usage

The pipeline is run using the `capcruncher pipeline` command. Ensure that you have a configuration file and the fastq files to process are in the current working directory.

```bash
capcruncher pipeline --cores <NUMBER OF CORES TO USE>
```

### HPC Cluster Usage

The pipeline can also be run on HPC clusters using a number of different job schedulers. See [here](cluster_config.md) for a quick overview on how to configure the pipeline to run on HPC clusters.

For further information see both the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) and the [Snakemake profile documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).

#### SLURM

This is a quick example of how to run the pipeline on a SLURM cluster. This is not a complete guide and you will need to modify the configuration to suit your cluster.

```bash
capcruncher pipeline -c <NUMBER OF CORES e.g. 20> --profile <NAME OF PROFILE OR PATH TO PROFILE>
```


### Singularity Usage (Recommended)

Containers have the advantage of their contents being fixed at the time of creation. This means that the pipeline will always run with the same versions of the software and aids reliablity and reproducibility. The pipeline can be run using singularity containers. This is the recommended method of running the pipeline.

The pipeline can be run using singularity containers using the `--use-singularity` option.

```bash
# Local mode
capcruncher pipeline --use-singularity --cores <NUMBER OF CORES TO USE>

# Cluster mode
capcruncher pipeline --use-singularity --cores <NUMBER OF CORES TO USE> --profile <NAME OF PROFILE OR PATH TO PROFILE>
```
