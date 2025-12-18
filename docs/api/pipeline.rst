API Reference: Pipeline
=======================

.. currentmodule:: qp_planck.qp_pipeline

High-level pipeline functions for running the full QuickPol workflow.

Main Functions
--------------

.. autofunction:: run_qp_pipeline

.. autofunction:: run_qp_from_yaml

Helper Functions
----------------

.. autofunction:: _load_yaml_config

.. autofunction:: _load_rimo_from_config

Examples
--------

Basic Pipeline Execution
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck import run_qp_pipeline
   
   run_qp_pipeline(
       detpairs=[("143A", "143A"), ("143A", "143B")],
       rimo_hfi="/data/RIMO_HFI.fits",
       blmfile="/data/beams/blm_{}.fits",
       outdir="./output",
       smax=6,
       release="npipe6v20"
   )

YAML Configuration
^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck import run_qp_from_yaml
   
   # Run from YAML file
   run_qp_from_yaml("config.yaml")
   
   # Or load config and override specific parameters
   run_qp_pipeline(
       detpairs=[],  # read from config
       config="config.yaml",
       overwrite=True  # override YAML value
   )

MPI Execution
^^^^^^^^^^^^^

.. code-block:: python

   # mpi_run.py
   from qp_planck import run_qp_pipeline
   
   # Define many detector pairs
   pairs = [
       ("100GHz", "100GHz"),
       ("143GHz", "143GHz"),
       ("217GHz", "217GHz"),
       # ... many more ...
   ]
   
   run_qp_pipeline(
       detpairs=pairs,
       rimo_hfi="/data/RIMO_HFI.fits",
       blmfile="/data/beams/blm_{}.fits",
       outdir="/scratch/beams"
   )

Execute with::

   mpirun -n 16 python mpi_run.py

Each MPI rank processes a subset of detector pairs automatically.

Custom RIMO Loading
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck.qp_pipeline import _load_rimo_from_config
   
   # Load from config dictionary
   config = {
       "rimo_lfi": "/data/RIMO_LFI.fits",
       "rimo_hfi": "/data/RIMO_HFI.fits"
   }
   
   RIMO = _load_rimo_from_config(config)
   print(f"Loaded {len(RIMO)} detectors")

Parameter Validation
^^^^^^^^^^^^^^^^^^^^

The pipeline validates inputs and provides helpful error messages:

.. code-block:: python

   from qp_planck import run_qp_pipeline
   
   try:
       run_qp_pipeline(
           detpairs=[("INVALID", "DETECTOR")],
           rimo_hfi="/data/RIMO_HFI.fits",
           blmfile="",
           outdir="./output"
       )
   except ValueError as e:
       print(f"Configuration error: {e}")

Common Patterns
^^^^^^^^^^^^^^^

**Test run before production:**

.. code-block:: python

   # Quick test with sparse sampling
   run_qp_pipeline(
       detpairs=[("143A", "143A")],
       rimo_hfi="/data/RIMO_HFI.fits",
       blmfile="/data/beams/blm_{}.fits",
       outdir="./test",
       test=True,  # Fast sparse sampling
       overwrite=True
   )
   
   # Full production run
   run_qp_pipeline(
       detpairs=[("143A", "143A")],
       rimo_hfi="/data/RIMO_HFI.fits",
       blmfile="/data/beams/blm_{}.fits",
       outdir="./production",
       test=False,  # Full sampling
       conserve_memory=True
   )

**Incremental processing:**

.. code-block:: python

   # Process different detector sets separately
   for freq in [100, 143, 217, 353]:
       pairs = [(f"{freq}GHz", f"{freq}GHz")]
       
       run_qp_pipeline(
           detpairs=pairs,
           rimo_hfi="/data/RIMO_HFI.fits",
           blmfile="/data/beams/blm_{}.fits",
           outdir=f"./beams_{freq}",
           release=f"npipe6v20_{freq}",
           overwrite=False  # Skip if already done
       )
