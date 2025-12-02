Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/>`_.

[Unreleased]
------------

Added
^^^^^
* Comprehensive Read the Docs style documentation
* API reference with examples for all modules
* Physics background section explaining QuickPol formalism
* YAML configuration guide
* Multiple usage examples
* Contributing guidelines

[0.1.0] - 2025
--------------

Added
^^^^^
* Initial release of qp-planck package
* Adapted QuickPol driver (``qp_hmap2mat.py``) from Planck NPIPE pipeline
* Beam matrix to FITS converter (``qp_mat2fits.py``)
* High-level pipeline wrapper (``qp_pipeline.py``) with YAML support
* Utility functions for RIMO loading and detector management
* Support for Planck LFI and HFI detectors
* MPI parallelization via mpi4py (optional)
* Gaussian beam fallback when beam multipoles unavailable
* Example YAML configuration file
* Helper script for generating detector pairs
* MIT License
* README with basic usage instructions
* GitHub Actions CI workflow for code style checks
* Read the Docs configuration
* Code coverage via codecov

Features
^^^^^^^^
* Compute effective beam window functions (Bℓ, Wℓ) for Planck detectors
* Full support for polarization (T, E, B coupling matrices)
* Asymmetric beam handling via beam multipoles
* Detector cross-talk and mismatch modeling
* Scanning strategy effects through hit maps
* Memory-efficient processing modes
* Python 3.9+ compatibility
* Type-hinted codebase
* Google-style docstrings

Technical Details
^^^^^^^^^^^^^^^^^
* Core dependencies: numpy, scipy, astropy, healpy, matplotlib
* Optional dependencies: PyYAML (for YAML configs), mpi4py (for parallel execution)
* Development tools: pytest, ruff, pre-commit, sphinx
* Package build: hatchling backend
* Testing: pytest with coverage reporting
* Documentation: Sphinx with RTD theme

Known Limitations
^^^^^^^^^^^^^^^^^
* Hard-coded paths for hit map directories (requires manual editing)
* Planet deconvolution is experimental
* Masking functionality is currently a stub
* Some functions have minimal documentation in source

[0.0.1] - 2024 (Pre-release)
-----------------------------

Changed
^^^^^^^
* Removed TOAST dependencies from original NPIPE code
* Restructured as standalone Python package
* Modernized packaging with pyproject.toml

---

**Note**: Version 0.1.0 is the first public release. Earlier development
was based on internal adaptations of the Planck NPIPE pipeline code.

See the `GitHub releases page <https://github.com/paganol/qp_planck/releases>`_
for more details on each version.
