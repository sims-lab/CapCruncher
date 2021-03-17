.. capture-c documentation master file, created by
   sphinx-quickstart on Wed Jul  1 11:37:07 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CCanalyser documentation
########################

CCanalyser is a package explicitly designed for processing Capture-C, Tri-C and Tiled-C data.
It consists of a configurable data processing pipeline and a supporting command line interface
to enable fine-grained control over the analysis. Unlike other pipelines that are designed to process
Hi-C or Capture-HiC data, the filtering steps in ccanalyser are specifically optimised Capture-C, Tri-C 
and Tiled-C. This package is effectively a python based re-implementation of `CCseqBasicS <https://github.com/Hughes-Genome-Group/CCseqBasicS>`_,
`Tri-C <https://github.com/oudelaar/TriC>`_, `Tiled-C <https://github.com/oudelaar/TiledC>`_ and 
`CaptureCompare <https://github.com/djdownes/CaptureCompare>`_ designed to provide an end-to-end
solution for processing data from these three assays.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installing.rst
   running_the_pipeline.rst
   faq_pipeline.rst
   cli.rst
   ccanalyser_module_documentation/modules.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
