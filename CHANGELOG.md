# Changelog

All notable changes to this project will be documented in this file.

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

## [0.2.2] - 2022-03-03

### Bug Fixes

- Fixes error when generating a report for a single sample([#149](https://github.com/sims-lab/CapCruncher/pull/149))
- Prevent cluster errors due to lack of allocated memory. ([#153](https://github.com/sims-lab/CapCruncher/pull/153))
- Enabled max memory adjustment for alignment filtering ([#156](https://github.com/sims-lab/CapCruncher/pull/156))
- Fix dask client shutdown related errors ([#154](https://github.com/sims-lab/CapCruncher/pull/154))
- Fix issue with missing report text file ([#157](https://github.com/sims-lab/CapCruncher/pull/157))

### Features

- Added CLI utility ti get viewpoint coordinates from capture oligos ([#144](https://github.com/sims-lab/CapCruncher/pull/144))
- Enables split directly to gzip compressed fastq to reduce the disk space required ([#146](https://github.com/sims-lab/CapCruncher/pull/146))
- Improve utility of the reporters-compare CLI module ([#145](https://github.com/sims-lab/CapCruncher/pull/145))
- Reduced disk space required by fastq files ([#150](https://github.com/sims-lab/CapCruncher/pull/150))
- Enabled biasing of duplicate removal towards trans reporters ([#151](https://github.com/sims-lab/CapCruncher/pull/151))

### Miscellaneous Tasks

- Cache conda env to speed up workflow ([#143](https://github.com/sims-lab/CapCruncher/pull/143))
- Publish development builds (test-pypi) and releases (pypi)  ([#152](https://github.com/sims-lab/CapCruncher/pull/152))
- Pre-release checks ([#159](https://github.com/sims-lab/CapCruncher/pull/159))
- Prepare for release 0.2.2

### Build

- Add issue and feature request templates ([#147](https://github.com/sims-lab/CapCruncher/pull/147))

## [0.2.0] - 2022-02-08

### Bug Fixes

- Fixed packaging long description. ([#115](https://github.com/sims-lab/CapCruncher/pull/115))
- Added missing dependencies (seaborn and trackhub) to setup.cfg
- Fixed bug after updating pyyaml to latest version ([#122](https://github.com/sims-lab/CapCruncher/pull/122))
- Re-partition reporter slices after filtering ([#124](https://github.com/sims-lab/CapCruncher/pull/124))
- Fixed help option not being displayed when running capcruncher pipeline ([#129](https://github.com/sims-lab/CapCruncher/pull/129))
- Fix issue with tasks going over their allotted number of cores. ([#133](https://github.com/sims-lab/CapCruncher/pull/133))
- Fixes error during deduplication when using gzip compression ([#134](https://github.com/sims-lab/CapCruncher/pull/134))
- Prevents excessive memory usage during duplicate removal ([#136](https://github.com/sims-lab/CapCruncher/pull/136))
- Fix link common cooler tables ([#137](https://github.com/sims-lab/CapCruncher/pull/137))
- Fixed an issue when no data exists for a viewpoint for a given sample  ([#139](https://github.com/sims-lab/CapCruncher/pull/139))

### Co-authored-by

- Kevin Rue-Albrecht <kevinrue67@gmail.com>
- Alsmith151 <49727900+alsmith151@users.noreply.github.com>
- Alsmith151 <49727900+alsmith151@users.noreply.github.com>

### Documentation

- Update README ([#120](https://github.com/sims-lab/CapCruncher/pull/120))

### Features

- Moved all configuration from setup.py to setup.cfg. ([#114](https://github.com/sims-lab/CapCruncher/pull/114))
- Reporter counting is now performed in parallel on separate partitions before collating. ([#117](https://github.com/sims-lab/CapCruncher/pull/117))
- Enables pileup normalisation using a set of regions supplied as a bed file ([#121](https://github.com/sims-lab/CapCruncher/pull/121))
- Enabled the use of custom filtering orders ([#119](https://github.com/sims-lab/CapCruncher/pull/119))
- Capability to normalise pileups (bedgraphs/bigwigs) by a set of supplied regions. ([#125](https://github.com/sims-lab/CapCruncher/pull/125))
- Expand the number of viewpoints that can be processed  ([#128](https://github.com/sims-lab/CapCruncher/pull/128))
- Enable optional compression during fastq split and deduplicate ([#131](https://github.com/sims-lab/CapCruncher/pull/131))
- Reduces disk space required by pipeline by removing intermediate files ([#135](https://github.com/sims-lab/CapCruncher/pull/135))
- Reduce disk space taken up by reporters (slices and counts) ([#138](https://github.com/sims-lab/CapCruncher/pull/138))
- Improved report generation and fixed multiple report related bugs ([#132](https://github.com/sims-lab/CapCruncher/pull/132))
- Revert without_cluster for reporter comparisons ([#140](https://github.com/sims-lab/CapCruncher/pull/140))

### REMOVED

- Get_yaml.py and test_run_params.py
- Validate_annotations.py