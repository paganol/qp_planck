API Reference: Utilities
=========================

.. currentmodule:: qp_planck.utilities

This module provides utility functions for working with Planck detector data,
including RIMO loading, detector lists, and filename helpers.

RIMO Loading
------------

.. autofunction:: load_RIMO

.. autoclass:: DetectorData
   :members:

Detector Lists
--------------

.. autofunction:: list_planck

Detector Weights
----------------

.. autodata:: detector_weights
   :annotation: = {detector: weight}
   
   Dictionary mapping detector names (strings) to white-noise weights (floats).
   Used for optimal map combination.
   
   Example::
   
      from qp_planck import detector_weights
      
      weight_143_1a = detector_weights["143-1"]
      print(f"Weight for 143-1 horn: {weight_143_1a:.2e}")

Filename Helpers
----------------

.. autofunction:: qp_file

Quaternion Operations
---------------------

.. autofunction:: quat_mult

Helper Classes
--------------

.. autoclass:: Timer
   :members:

Constants
---------

.. autodata:: degree
   :annotation: = π/180
   
   Conversion factor from degrees to radians.

.. autodata:: SPINROT
   :annotation: = [0, 0, 0, 1]
   
   Bore-sight rotation quaternion. Set to identity (no rotation) by default.

Examples
--------

Loading RIMO and Accessing Detector Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck import load_RIMO
   
   # Load HFI RIMO
   RIMO = load_RIMO("RIMO_HFI_npipe5v16_symmetrized.fits")
   
   # Access detector properties
   det = RIMO["143-1a"]
   print(f"Detector: {det.name}")
   print(f"FWHM: {det.fwhm:.2f} arcmin")
   print(f"Polarization angle: {det.psi_pol:.2f} degrees")
   print(f"Polarization efficiency ε: {det.epsilon:.4f}")
   print(f"NET: {det.net:.2f} μK√s")
   
   # Compute ρ = (1-ε)/(1+ε)
   rho = (1 - det.epsilon) / (1 + det.epsilon)
   print(f"ρ: {rho:.4f}")

Listing Detectors
^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck import list_planck
   
   # All 143 GHz detectors
   dets_143 = list_planck("143GHz")
   print(f"143 GHz: {dets_143}")
   
   # Subset A
   dets_143A = list_planck("143A")
   print(f"143A: {dets_143A}")
   
   # Single horn (combines PSBs)
   horn_143_1 = list_planck("143-1")
   print(f"143-1 horn: {horn_143_1}")
   
   # All Planck detectors
   all_planck = list_planck("Planck")
   print(f"Total Planck detectors: {len(all_planck)}")

Generating Filenames
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck import qp_file
   
   filename = qp_file(
       outdir="./beam_output",
       dets=("143A", "143B"),
       lmax=4096,
       smax=6,
       angle_shift=0.0,
       full=True,
       pconv="cmbfast",
       release="npipe6v20",
       rhobeam="IMO",
       rhohit="IMO"
   )
   print(filename)
   # Output: ./beam_output/beam_matrix_143Ax143B_l4096_s6_A000_cmbfast_1_rbIMO_rhIMO.npz

Working with MPI
^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck import load_RIMO
   
   # RIMO loading with MPI (automatic broadcast)
   try:
       from mpi4py import MPI
       comm = MPI.COMM_WORLD
   except ImportError:
       comm = None
   
   RIMO = load_RIMO("RIMO_HFI.fits", comm=comm)
   
   # All ranks now have the RIMO
   if comm is None or comm.rank == 0:
       print(f"Loaded {len(RIMO)} detectors")
