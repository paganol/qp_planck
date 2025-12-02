Installation
============

Requirements
------------

Python Version
^^^^^^^^^^^^^^

**qp-planck** requires Python 3.9 or later.

Dependencies
^^^^^^^^^^^^

Core dependencies (automatically installed):

* **numpy** ≥ 1.20 — Array operations and numerical computing
* **scipy** ≥ 1.6 — Interpolation and optimization
* **astropy** ≥ 5.0 — FITS file I/O and astronomical utilities
* **healpy** ≥ 1.15 — HEALPix map operations and spherical harmonics
* **matplotlib** ≥ 3.9.4 — Plotting and visualization

Recommended dependencies:

* **PyYAML** — YAML configuration file parsing (required for ``run_qp_from_yaml``)

Optional dependencies:

* **mpi4py** ≥ 3.0 — MPI parallelization (recommended for large-scale runs)

Development dependencies (for contributing):

* **pytest** ≥ 8.4.2 — Testing framework
* **pytest-cov** ≥ 7.0.0 — Code coverage
* **ruff** ≥ 0.14.0 — Linting and formatting
* **pre-commit** ≥ 4.3.0 — Git hooks for code quality
* **sphinx** ≥ 7.0.0 — Documentation generation
* **sphinx-rtd-theme** ≥ 2.0.0 — Read the Docs theme
* **sphinx-autodoc-typehints** ≥ 2.0.0 — Type hints in documentation

Installation Methods
--------------------

From Source
^^^^^^^^^^^

Clone the repository and install in editable mode:

.. code-block:: bash

   git clone https://github.com/paganol/qp_planck.git
   cd qp_planck
   
   # Using uv (recommended)
   uv sync --all-extras
   
   # Or using pip
   pip install -e .

This allows you to modify the code and see changes immediately without reinstalling.

With MPI Support
^^^^^^^^^^^^^^^^

To enable MPI parallelization:

.. code-block:: bash

   pip install -e ".[mpi]"

This installs **mpi4py** in addition to the core dependencies.

For Developers
^^^^^^^^^^^^^^

Install with all development tools:

.. code-block:: bash

   pip install -e ".[dev]"

Or install everything (MPI + development):

.. code-block:: bash

   pip install -e ".[mpi,dev]"

Using uv (Recommended)
^^^^^^^^^^^^^^^^^^^^^^

If you have `uv <https://github.com/astral-sh/uv>`_ installed (recommended):

.. code-block:: bash

   # Install all dependencies including optional ones
   uv sync --all-extras
   
   # Or install specific extras
   uv sync --extra dev     # Development tools
   uv sync --extra mpi     # MPI support

This respects the ``uv.lock`` file for reproducible builds and is the preferred 
installation method for development.

Verification
------------

Test your installation:

.. code-block:: python

   import qp_planck
   from qp_planck import load_RIMO, list_planck
   
   # List all Planck detectors
   print(list_planck("Planck"))
   
   # Check if MPI is available
   try:
       from mpi4py import MPI
       print(f"MPI available with {MPI.COMM_WORLD.size} processes")
   except ImportError:
       print("MPI not available (serial mode only)")

System Requirements
-------------------

**Disk Space**

* Package installation: ~10 MB
* Input data (RIMO, beams, hitmaps): varies, typically 1-10 GB
* Output beam matrices: ~100 MB - 1 GB depending on detector coverage

**Memory**

* Minimum: 4 GB RAM
* Recommended: 16+ GB RAM for full-frequency beam matrix calculations
* Use ``conserve_memory=True`` option for memory-constrained systems

**Computing**

* Single-core execution: possible but slow for full detector sets
* MPI parallelization: highly recommended for production runs
* Typical HPC usage: 32-128 cores for full Planck coverage

Troubleshooting
---------------

HEALPix Installation Issues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If **healpy** fails to install, ensure you have a working C compiler:

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt-get install build-essential
   
   # macOS
   xcode-select --install

MPI Compilation Problems
^^^^^^^^^^^^^^^^^^^^^^^^^

For **mpi4py**, you need an MPI implementation:

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt-get install libopenmpi-dev
   
   # macOS (with Homebrew)
   brew install open-mpi

Then install mpi4py:

.. code-block:: bash

   pip install mpi4py

Import Errors
^^^^^^^^^^^^^

If you get ``ModuleNotFoundError`` after installation:

1. Check you installed in the correct Python environment
2. Verify with: ``pip list | grep qp-planck``
3. Try reinstalling: ``pip install --force-reinstall -e .``

Next Steps
----------

* :doc:`quickstart` — Run your first beam window calculation
* :doc:`yaml_config` — Configure the pipeline with YAML
* :doc:`examples` — Explore example workflows
