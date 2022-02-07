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

###  v0.2.0 - (2022-02-07)

### Bug Fixes

* **CLI:** Fixed help option not being displayed when running capcruncher pipeline ([#129](https://github.com/sims-lab/CapCruncher/issues/129)) ([9f09093](https://github.com/sims-lab/CapCruncher/commit/9f090935f3c20c5d78e01ab4f5b0248b325ee341))
* **alignment_deduplication:** Prevents excessive memory usage during duplicate removal ([#136](https://github.com/sims-lab/CapCruncher/issues/136)) ([b175978](https://github.com/sims-lab/CapCruncher/commit/b17597884164ed074782370637a81732390ac48c))
* **packaging:** Fixed bug after updating pyyaml to latest version ([#122](https://github.com/sims-lab/CapCruncher/issues/122)) ([7d76b5f](https://github.com/sims-lab/CapCruncher/commit/7d76b5f4976fe3c6f1bc09989df3db28c12ecce3))
* Added missing dependencies (seaborn and trackhub) to setup.cfg ([550a882](https://github.com/sims-lab/CapCruncher/commit/550a882af5e131c04b5d45bf0430ecc50ce15310))
* Fixed packaging long description. ([#115](https://github.com/sims-lab/CapCruncher/issues/115)) ([6f716d1](https://github.com/sims-lab/CapCruncher/commit/6f716d182de705146333206a38d9c791de1a9227))
* **pipeline:** 
* Fixes issue with tasks going over their allotted number of cores. ([#133](https://github.com/sims-lab/CapCruncher/issues/133)) ([27cd193](https://github.com/sims-lab/CapCruncher/commit/27cd193c207409b96a0b28c079b9d689daaa61ee))
* Fixes error during deduplication when using gzip compression ([#134](https://github.com/sims-lab/CapCruncher/issues/134)) ([01ff56b](https://github.com/sims-lab/CapCruncher/commit/01ff56b88558af486d11b9f7544c8c5c6ca9f002))
* Re-partition reporter slices after filtering ([#124](https://github.com/sims-lab/CapCruncher/issues/124)) ([db72c56](https://github.com/sims-lab/CapCruncher/commit/db72c56875c13ed2762d44e916a8ed66f73324cc))
* **reporters_compare:** Fixed an issue when no data exists for a viewpoint for a given sample  ([#139](https://github.com/sims-lab/CapCruncher/issues/139)) ([e720029](https://github.com/sims-lab/CapCruncher/commit/e7200299bf2453e719f28f95ed3658e9570b7ad5))
* **storage:** Fix link common cooler tables ([#137](https://github.com/sims-lab/CapCruncher/issues/137)) ([4836fbe](https://github.com/sims-lab/CapCruncher/commit/4836fbe8e46ad268dda6d05f27104789f0c46e0d))


### Features

* **CLI:** Enables pileup normalisation using a set of regions supplied as a bed file ([#121](https://github.com/sims-lab/CapCruncher/issues/121)) ([9c587ff](https://github.com/sims-lab/CapCruncher/commit/9c587ff1a60f009c0b990952361810d61376a1c7))
* **pipeline**: Expanded the number of viewpoints that can be processed  ([#128](https://github.com/sims-lab/CapCruncher/issues/128)) ([8fcb576](https://github.com/sims-lab/CapCruncher/commit/8fcb57657f108d78cdbb1e255a5eb85b7cb3e860))
* **packaging:**  Moved all configuration from setup.py to setup.cfg. ([#114](https://github.com/sims-lab/CapCruncher/issues/114)) ([4835da4](https://github.com/sims-lab/CapCruncher/commit/4835da44157132feda38e299bf9c67ca297c3d2d))
* **pipeline:** 
* Capability to normalise pileups (bedgraphs/bigwigs) by a set of supplied regions. ([#125](https://github.com/sims-lab/CapCruncher/issues/125)) ([bab07ea](https://github.com/sims-lab/CapCruncher/commit/bab07eac1e524020d24c745dd88b749173d9d440))
* Enable optional compression during fastq split and deduplicate ([#131](https://github.com/sims-lab/CapCruncher/issues/131)) ([0c32b73](https://github.com/sims-lab/CapCruncher/commit/0c32b7320fcff5d95145a406996e9baf9f7aeebd))
* Enabled the use of custom filtering orders ([#119](https://github.com/sims-lab/CapCruncher/issues/119)) ([b57ebe8](https://github.com/sims-lab/CapCruncher/commit/b57ebe886fc767b8dcb12c7dfc45dd2e9a1ea1b3))
* Reduces disk space required by pipeline by removing intermediate files ([#135](https://github.com/sims-lab/CapCruncher/issues/135)) ([d6c4302](https://github.com/sims-lab/CapCruncher/commit/d6c4302a27c14b965c531b11242ef6dd152fc1a1))
* Reporter counting is now performed in parallel on separate partitions before collating. ([#117](https://github.com/sims-lab/CapCruncher/issues/117)) ([aae5356](https://github.com/sims-lab/CapCruncher/commit/aae5356d6268e71ae777ffb31fcbd98e76ccd8c2))
* Revert without_cluster for reporter comparisons ([#140](https://github.com/sims-lab/CapCruncher/issues/140)) ([f847d28](https://github.com/sims-lab/CapCruncher/commit/f847d282f556d336be2a66023aced8c8dd082551))
* **storage:** Reduce disk space taken up by reporters (slices and counts) ([#138](https://github.com/sims-lab/CapCruncher/issues/138)) ([7659a8c](https://github.com/sims-lab/CapCruncher/commit/7659a8c3fee15ec94c107313d16ce9c831f4ffbf))



