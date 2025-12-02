Quick Start Guide
=================

This guide will walk you through computing beam window functions for Planck
detectors using **qp-planck**.

Prerequisites
-------------

Before starting, ensure you have:

1. **Installed qp-planck** (see :doc:`installation`)

   .. code-block:: bash
   
      # Using uv (recommended)
      uv sync --all-extras
      
      # Or using pip
      pip install -e .

2. **Input data files**:

   * RIMO files: ``RIMO_LFI_*.fits`` and ``RIMO_HFI_*.fits``
   * Beam multipoles: ``blm_<detector>.fits`` for each detector
   * Hit maps and spin moments: ``polmoments_<detector>.fits`` files

These files are typically from the Planck NPIPE data release. Contact the
Planck team or check the `Planck Legacy Archive <https://pla.esac.esa.int/>`_
for access.

Basic Workflow
--------------

The **qp-planck** pipeline has three main steps:

1. **Load RIMO** — Read detector properties (angles, FWHM, noise, etc.)
2. **Compute beam matrices** — Combine beams + hitmaps → NPZ files
3. **Convert to FITS** — Extract :math:`B_\ell` and :math:`W_\ell` window functions

Step 1: Quick Test Run
-----------------------

Start with a simple Python script to verify your setup:

.. code-block:: python

   from qp_planck import load_RIMO, list_planck
   
   # Load Planck detector properties
   RIMO = load_RIMO("path/to/RIMO_HFI_npipe5v16_symmetrized.fits")
   
   # List detectors in a frequency channel
   detectors_143 = list_planck("143GHz")
   print(f"143 GHz detectors: {detectors_143}")
   
   # Access properties of a detector
   det = "143-1a"
   print(f"Detector {det}:")
   print(f"  FWHM: {RIMO[det].fwhm:.2f} arcmin")
   print(f"  Polarization efficiency: {RIMO[det].epsilon:.3f}")
   print(f"  Noise (NET): {RIMO[det].net:.2f} μK√s")

Expected output::

   143 GHz detectors: ['143-1a', '143-1b', '143-2a', '143-2b', ..., '143-7']
   Detector 143-1a:
     FWHM: 7.22 arcmin
     Polarization efficiency: 0.002
     Noise (NET): 40.12 μK√s

Step 2: Compute Beam Matrices
------------------------------

Now compute the beam matrix for a single detector pair:

.. code-block:: python

   from qp_planck import run_qp_pipeline
   
   run_qp_pipeline(
       detpairs=[("143-1a", "143-1a")],  # Auto-spectrum for detector 143-1a
       rimo_hfi="path/to/RIMO_HFI_npipe5v16_symmetrized.fits",
       blmfile="path/to/beams/blm_{}.fits",
       outdir="./beam_output",
       smax=6,
       spin_ref="Pxx",
       blm_ref="Dxx",
       release="npipe6v20",
       test=False,  # False = full sampling
       overwrite=True
   )

This will:

* Load the RIMO for detector ``143-1a``
* Read beam multipoles from ``blm_143-1a.fits``
* Load hit maps/spin moments
* Compute the beam matrix
* Save results to ``beam_output/beam_matrix_143-1ax143-1a_*.npz``
* Convert to FITS: ``Bl_npipe6v20_143-1ax143-1a.fits`` and ``Wl_*.fits``

.. note::
   The first run may take several minutes depending on your system.
   Set ``test=True`` for faster debugging (sparse pixel sampling).

Step 3: Inspect Results
-----------------------

Load and plot the beam window function:

.. code-block:: python

   import matplotlib.pyplot as plt
   from astropy.io import fits
   import numpy as np
   
   # Load the temperature beam window
   with fits.open("beam_output/Bl_npipe6v20_143-1ax143-1a.fits") as hdul:
       bl = hdul[1].data['TEMPERATURE']
   
   # Plot
   ell = np.arange(len(bl))
   plt.figure(figsize=(10, 6))
   plt.semilogy(ell, bl, label='143-1a beam')
   plt.xlabel(r'Multipole $\ell$')
   plt.ylabel(r'$B_\ell$')
   plt.title('Planck 143 GHz Beam Window Function')
   plt.grid(True, alpha=0.3)
   plt.legend()
   plt.xlim(0, 3000)
   plt.ylim(1e-8, 1)
   plt.savefig('beam_window_143-1a.png', dpi=150)
   plt.show()

You should see a smooth curve peaking at :math:`B_0 = 1` and decreasing
with :math:`\ell`, reflecting beam damping of small-scale power.

Running with YAML Configuration
--------------------------------

For production runs with many detector pairs, use a YAML config file:

.. code-block:: yaml

   # beam_config.yaml
   
   # RIMO files
   rimo_hfi: /path/to/RIMO_HFI_npipe5v16_symmetrized.fits
   
   # Beam multipoles
   blmfile: /path/to/beams/blm_{}.fits
   
   # Output directory
   outdir: ./beam_output
   
   # Physics parameters
   smax: 6
   spin_ref: Pxx
   blm_ref: Dxx
   release: npipe6v20
   
   # Computational settings
   test: false
   conserve_memory: true
   overwrite: false
   
   # FITS output options
   blfile: true      # Scalar B_l(T)
   blTEBfile: true   # B_l for T, E, B
   wlfile: true      # Full W_l matrices
   
   # Detector pairs (auto-spectra for 143 GHz)
   detpairs:
     - ['143-1a', '143-1a']
     - ['143-1b', '143-1b']
     - ['143-2a', '143-2a']
     - ['143-2b', '143-2b']
     - ['143-3a', '143-3a']
     - ['143-3b', '143-3b']
     - ['143-4a', '143-4a']
     - ['143-4b', '143-4b']

Run it with:

.. code-block:: bash

   python -m qp_planck.qp_pipeline beam_config.yaml

Or from Python:

.. code-block:: python

   from qp_planck import run_qp_from_yaml
   
   run_qp_from_yaml("beam_config.yaml")

MPI Parallel Execution
----------------------

For large-scale runs, use MPI to distribute detector pairs across cores:

.. code-block:: bash

   mpirun -n 8 python -m qp_planck.qp_pipeline beam_config.yaml

This will:

* Split the detector pairs among 8 MPI ranks
* Each rank processes its subset independently
* Results are written to separate files (no communication needed)

.. tip::
   Choose ``n`` ≤ number of detector pairs for best efficiency.
   Overhead is minimal, so even modest parallelism (4-8 cores) helps significantly.

Understanding the Outputs
--------------------------

After a successful run, you'll have:

**NPZ Beam Matrices** (``beam_matrix_*.npz``):

Intermediate files containing:

* ``beam_mat``: dictionary of beam matrices for TT, EE, BB, TE, TB, EB
* ``hit_mat``: hit matrix for all spin components
* ``lmax``, ``smax``, detector info, metadata

These are typically not needed for analysis but can be useful for debugging.

**FITS Beam Windows** (``Bl_*.fits``):

Scalar beam window :math:`B_\ell` for temperature:

* Extension 1: column ``TEMPERATURE`` with :math:`B_\ell` for :math:`\ell = 0, \ldots, \ell_{\max}`

Compatible with HEALPix tools like ``synfast`` and PolSpice.

**FITS TEB Beam Windows** (``Bl_TEB_*.fits``):

Separate :math:`B_\ell` for T, E, B:

* Extension 1: columns ``T``, ``E``, ``B``

Used for approximate diagonal beam correction in polarization analyses.

**FITS Window Matrices** (``Wl_*.fits``):

Full coupling matrices :math:`W_\ell`:

* Extension 1 (TT): 9 columns for :math:`W_\ell^{TT \to XY}`
* Extension 2 (EE): 9 columns for :math:`W_\ell^{EE \to XY}`
* Extension 3 (BB): 9 columns for :math:`W_\ell^{BB \to XY}`
* Extension 4 (TE): 9 columns for :math:`W_\ell^{TE \to XY}`

These account for full T/E/B coupling and are required for high-precision
power spectrum estimation.

Common Issues
-------------

**Missing Input Files**

If you see errors like::

   FileNotFoundError: polmoments_143-1a.fits not found

Ensure the hit maps/spin moments are in the directory specified by the
code's internal ``hitgrpfull`` variable. You may need to edit
``qp_planck/qp_hmap2mat.py`` to point to your data location.

**Memory Errors**

For large detector sets, use:

.. code-block:: python

   run_qp_pipeline(..., conserve_memory=True, test=True)

This reduces memory usage and uses sparse pixel sampling for testing.

**Slow Execution**

* Use MPI: ``mpirun -n <cores> python -m qp_planck.qp_pipeline config.yaml``
* Reduce :math:`\ell_{\max}`: not recommended, but possible by editing detector NSIDE
* Skip cross-spectra: only compute auto-spectra (diagonal terms)

Next Steps
----------

* :doc:`yaml_config` — Full YAML configuration options
* :doc:`physics_background` — Understand the QuickPol formalism
* :doc:`examples` — Advanced use cases (cross-spectra, frequency maps, etc.)
* :doc:`api/pipeline` — API reference for ``run_qp_pipeline``
