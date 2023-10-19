# CapCruncher Tips and Tricks

## Interruptions to the pipeline (e.g. system failure, manual interruption)

### Restarting the pipeline after an interruption

If the pipeline is interrupted, it can be restarted by simply running the pipeline command again e.g. `capcruncher pipeline -c 1`.

CapCruncher will detect which steps have already been completed and will skip them. This is useful if the pipeline is interrupted due to a system failure or if you want to add more samples to the pipeline.

### Unlocking the working directory

Snakemake locks the working directory during the pipeline run. If the pipeline is interrupted, the working directory will remain locked and will not restart. To unlock the working directory, run:

``` bash
capcruncher pipeline --unlock
```

## Interuptions to the pipeline (e.g. error in pipeline)

Pipeline errors very frequently are found in a few major areas:

### No fastq files found

The pipeline cannot find the fastq files to process (e.g. the files are not in the current working directory or are not named correctly) this will cause an error like this:

``` bash
2023-08-03 11:56:17.857 | ERROR    | capcruncher.pipeline.utils:from_files:178 - No fastq files found.
ValueError in file /ceph/home/a/asmith/software/CapCruncher/capcruncher/pipeline/workflow/Snakefile, line 30:
No fastq files found.
  File "/ceph/home/a/asmith/software/CapCruncher/capcruncher/pipeline/workflow/Snakefile", line 30, in <module>
  File "/ceph/home/a/asmith/software/CapCruncher/capcruncher/pipeline/utils.py", line 179, in from_files
```

To correct this, ensure that the fastq files are in the current working directory and are named correctly. See the [pipeline guide](pipeline.md#setting-up-the-input-directory) for more information on how to name the fastq files.

### Indicies not found or specified incorrectly

If the pipeline cannot find the reference genome files this will cause an error at the alignment step. Please ensure that the reference genome (`aligner_indicies`) are added to the [config file](pipeline.md#pipeline-configuration).

```yaml
aligner_indicies: "/ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/chr14_bowtie2_indicies/bt2"
```

This refers to the bowtie2 index files here:

``` bash
tree "/ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/chr14_bowtie2_indicies/"
/ceph/home/a/asmith/software/CapCruncher/tests/data/data_for_pipeline_run/chr14_bowtie2_indicies/
|-- bt2.1.bt2
|-- bt2.2.bt2
|-- bt2.3.bt2
|-- bt2.4.bt2
|-- bt2.rev.1.bt2
|-- bt2.rev.2.bt2
```

### Incorrect viewpoint coordinates

Viewpoint coordinate errors fall into two categories:

#### Viewpoint names contain invalid characters

Including special characters e.g. "\/*?" in the viewpoint name will prevent the pipeline from running to completion. Viewpoint names should only contain alphanumeric characters (`A-Za-z0-9`), hyphens (`-`) and underscores (`_`). No other stipulations are placed on the viewpoint name.

#### Viewpoint coordinates are incorrect for the supplied reference genome

!!! warning
    Errors in the viewpoint coordinates can be difficult to spot as pipeline errors will occur further downstream. The initial error occurs at the filtering step but the pipeline will continue to run

    In the future the presence of valid viewpoints will be confirmed during the pipeline run but for now it is up to the user to ensure that the viewpoint coordinates are correct.


Viewpoint coordinates are supplied in the [config file](pipeline.md#pipeline-configuration) as a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file and should be checked against the reference genome.

##### Capture-C and Tri-C experiments

The viewpoint coordinates specified should be that of the restriction fragment containing the viewpoint. This is usually 1-4 kb in size but the coordinates do not have to be exact (e.g. +/- 100 bp is fine).

##### Tiled-C experiments

The viewpoint coordinates should contain all restriction fragments that have been captured by the tiled oligos. Again these do not have to be basepair level accurate but should be as close as possible.


















## Adding additional Snakemake options

Additional Snakemake options can be passed to the pipeline command by just adding them to the end of the command. For example, to run the pipeline with 8 cores and prevent the pipeline from removing intermediate files, run:

``` bash
capcruncher pipeline --cores 8 --notemp
```

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for a list of available options.
