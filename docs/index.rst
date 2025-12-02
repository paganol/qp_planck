qp-planck: QuickPol Beam Windows for Planck NPIPE
==================================================

Welcome to the **qp-planck** documentation! This package provides QuickPol-based
beam window function computation for the Planck NPIPE data release.

.. warning::
   **AI-Generated Documentation**: This documentation was automatically generated 
   using AI assistance. While comprehensive, it may contain inaccuracies or be 
   outdated. For the most accurate and up-to-date information, please refer to the 
   ``README.md`` and source code comments in the repository.

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT

.. image:: https://readthedocs.org/projects/qp-planck/badge/?version=latest
   :target: https://qp-planck.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

Overview
--------

The **qp-planck** package computes effective beam window functions for Planck
satellite detectors using the QuickPol formalism. It processes:

* **Reduced Instrument Model (RIMO)** files containing detector properties
* **Beam multipoles** (:math:`b_{\ell m}`) from optical modeling
* **Hit maps** and **spin moments** from scanning strategy
* Produces **beam matrices** and **window functions** (:math:`B_\ell`, :math:`W_\ell`)

The pipeline accounts for:

* Asymmetric beams and cross-polarization
* Detector orientation and polarization angle
* Scanning strategy effects
* Temperature-to-polarization leakage

.. warning::
   This code is **adapted from the Planck NPIPE pipeline**:
   
   https://github.com/planck-npipe/toast-npipe
   
   The original implementation relied on TOAST and internal Planck tooling.
   This version removes TOAST dependencies and runs as a standalone package.

Quick Start
-----------

Installation
^^^^^^^^^^^^

Install from the repository:

.. code-block:: bash

   git clone https://github.com/paganol/qp_planck.git
   cd qp_planck
   
   # Using uv (recommended)
   uv sync --all-extras
   
   # Or using pip
   pip install -e .

Required dependencies:

* numpy ≥ 1.20
* scipy ≥ 1.6
* astropy ≥ 5.0
* healpy ≥ 1.15
* matplotlib ≥ 3.9.4

Optional dependencies:

* PyYAML (for YAML configuration files)
* mpi4py ≥ 3.0 (for MPI parallelization)

Basic Usage
^^^^^^^^^^^

The easiest way to run the pipeline is with a YAML configuration file:

.. code-block:: bash

   python -m qp_planck.qp_pipeline scripts/example_planck.yaml

Or from Python:

.. code-block:: python

   from qp_planck import run_qp_pipeline

   run_qp_pipeline(
       detpairs=[("143A", "143A"), ("143A", "143B")],
       rimo_lfi="path/to/RIMO_LFI.fits",
       rimo_hfi="path/to/RIMO_HFI.fits",
       blmfile="path/to/beams/blm_{}.fits",
       outdir="output/",
       smax=6,
       release="npipe6v20"
   )

What's Inside
-------------

The package contains four main modules:

:mod:`qp_planck.utilities`
   RIMO loading, detector lists, weights, filename helpers

:mod:`qp_planck.qp_hmap2mat`
   QuickPol driver: builds beam matrices from hitmaps and beam multipoles

:mod:`qp_planck.qp_mat2fits`
   Converts beam matrices (NPZ) to FITS window functions

:mod:`qp_planck.qp_pipeline`
   High-level pipeline wrapper with YAML configuration support

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   physics_background
   configuration
   yaml_config
   examples

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/utilities
   api/hmap2mat
   api/mat2fits
   api/pipeline

.. toctree::
   :maxdepth: 1
   :caption: Additional Information

   references
   contributing
   changelog

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
