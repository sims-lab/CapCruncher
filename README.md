# CapCruncher

[![Documentation Status](https://readthedocs.org/projects/capcruncher/badge/?version=latest)](https://capcruncher.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/sims-lab/CapCruncher/branch/master/graph/badge.svg?token=RHIGNMGX09)](https://codecov.io/gh/sims-lab/CapCruncher)
![CI](https://github.com/sims-lab/capture-c/actions/workflows/python-template.yml/badge.svg)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/capcruncher/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![DOI](https://zenodo.org/badge/224631087.svg)](https://zenodo.org/badge/latestdoi/224631087)
[![Downloads](https://pepy.tech/badge/capcruncher)](https://pepy.tech/project/capcruncher)

## Analysis software for Capture-C, Tri-C and Tiled-C data.

CapCruncher is a tool designed to automate the processing of Capture-C, Tri-C and Tiled-C data from FASTQ files, the package is written in python and  consists of an end-to-end data processing pipline together with a supporting command line interface to enable finer grained control. The pipeline provided is fast, robust and  scales from a laptop to a computational cluster. 

For further information see the [documentation](https://capcruncher.readthedocs.io/en/latest/)

## Changelog

## [0.2.3] - 2022-08-19

### Bug Fixes

- Fixes Tiled-C filtering error due to typo in index specification for remove_dual_capture_fragments
- Fixed bug when annotating Tiled-C data ([#166](https://github.com/sims-lab/CapCruncher/pull/166)) with the sorted option that caused no data to be annotated as as a viewpoint or exclusion.
- Fixed a bug with Tiled-C slice filtering ([#166](https://github.com/sims-lab/CapCruncher/pull/166)) that caused slices to be erronously filtered.
- Fixed a bug with counting reporters in batches (this occurs when counting >1x106 slices per viewpoint) ([#166](https://github.com/sims-lab/CapCruncher/pull/166)) 
- Fixed a bug when merging and storing Tri-C or Tiled-C data ([#166](https://github.com/sims-lab/CapCruncher/pull/166)) using "capcruncher reporters store merge" or "capcruncher reporters store bins". This functionally has been re-written and now appears to work correctly.
- Fixed a bug with plotting matrices using capcrunchers plotting capabilities ([#166](https://github.com/sims-lab/CapCruncher/pull/166)). 
- Fixes bug where all slices are removed from parquet files (reporter) outputs due to an upgraded dask version (the pyarrow-dataset engine has been removed). This corrects the all dataframes are empty error occurring while generating reporter statistics.

### Features

- Added option to normalised CCMatrix data using ice followed by a scaling factor.
- Reporter merging (paraquet) has been re-written to use pyarrow directly and is now faster and better able to split datasets into smaller files for more efficient parallel querying.


### Miscellaneous Tasks

- Updated version to 0.2.3 ([#165](https://github.com/sims-lab/CapCruncher/pull/165))
- Pin conda dependency versions to speed up environment solver ([#167](https://github.com/sims-lab/CapCruncher/pull/167))