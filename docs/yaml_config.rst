YAML Configuration
==================

The **qp-planck** pipeline can be fully configured using YAML files, making it
easy to reproduce runs and manage complex detector combinations.

Configuration File Structure
-----------------------------

A complete YAML configuration file has the following sections:

.. code-block:: yaml

   # RIMO Files
   rimo_lfi: /path/to/RIMO_LFI.fits
   rimo_hfi: /path/to/RIMO_HFI.fits
   # OR use a single merged RIMO:
   # rimo: /path/to/RIMO_merged.fits
   
   # Beam Multipoles
   blmfile: /path/to/beams/blm_{}.fits
   
   # I/O Paths
   outdir: ./beam_output
   indir: ./beam_output  # where to read NPZ files for mat2fits
   
   # Physics Parameters
   smax: 6
   spin_ref: Pxx
   blm_ref: Dxx
   angle_shift: 0.0
   
   # Data Release
   release: npipe6v20
   rhobeam: IMO
   rhohit: IMO
   
   # Computational Settings
   test: false
   conserve_memory: true
   overwrite: false
   
   # Output Options
   blfile: true
   blTEBfile: true
   wlfile: true
   do_plot: false
   
   # Detector Pairs
   detpairs:
     - ['143A', '143A']
     - ['143A', '143B']
     - ['143B', '143B']

Configuration Parameters
------------------------

RIMO Files
^^^^^^^^^^

Specify paths to Reduced Instrument Model FITS files:

.. code-block:: yaml

   # Option 1: Separate LFI and HFI RIMOs (recommended)
   rimo_lfi: /path/to/RIMO_LFI_npipe5_symmetrized.fits
   rimo_hfi: /path/to/RIMO_HFI_npipe5v16_symmetrized.fits
   
   # Option 2: Single merged RIMO
   rimo: /path/to/RIMO_Planck_merged.fits

The RIMO files contain detector properties:

* Focal plane positions (``phi_uv``, ``theta_uv``, ``psi_uv``)
* Polarization angles (``psi_pol``)
* Polarization efficiency (``epsilon``)
* Beam FWHM (``fwhm``)
* Noise properties (``f_samp``, ``f_knee``, ``alpha``, ``net``)

Beam Multipoles
^^^^^^^^^^^^^^^

Template for beam multipole FITS file paths:

.. code-block:: yaml

   blmfile: /path/to/beams/blm_{}.fits

The ``{}`` placeholder is replaced by the detector name (e.g., ``143-1a``).

**Alternative: Gaussian beams**

Set ``blmfile: ""`` to use Gaussian beams based on RIMO FWHM instead of
reading measured :math:`b_{\ell m}`.

Input/Output Directories
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

   outdir: ./beam_output  # Where NPZ and FITS files are written
   indir: ./beam_output   # Where NPZ files are read from (for mat2fits step)

Typically ``outdir == indir`` unless you're reprocessing existing NPZ files.

Physics Parameters
^^^^^^^^^^^^^^^^^^

``smax``
   **Maximum spin moment** (integer, default: 6)
   
   Controls the expansion order in spin-weighted harmonics. Higher values
   increase accuracy but cost more memory/time. For Planck, ``smax=6`` is
   sufficient.

``spin_ref``
   **Spin moment reference frame** (``'Pxx'`` or ``'Dxx'``, default: ``'Pxx'``)
   
   * ``'Pxx'``: Detector polarization-sensitive direction
   * ``'Dxx'``: Detector :math:`uv` frame
   
   Must match the reference used when generating the hit maps / spin moments.

``blm_ref``
   **Beam multipole reference frame** (``'Pxx'`` or ``'Dxx'``, default: ``'Dxx'``)
   
   Reference frame in which the input :math:`b_{\ell m}` are defined.

``angle_shift``
   **Additional detector angle rotation** (float, degrees, default: 0.0)
   
   Applied to all detectors. Used for systematics studies.

Polarization Models
^^^^^^^^^^^^^^^^^^^

``rhobeam``
   **Beam polarization efficiency model** (``'IMO'`` or ``'Ideal'``, default: ``'IMO'``)
   
   * ``'IMO'``: Use measured :math:`\rho = (1-\epsilon)/(1+\epsilon)` from RIMO
   * ``'Ideal'``: Assume :math:`\rho = 1` for PSBs, :math:`\rho = 0` for SWBs

``rhohit``
   **Hit matrix polarization model** (``'IMO'`` or ``'Ideal'``, default: ``'IMO'``)
   
   Same options as ``rhobeam``, applied to the hit matrix calculation.

Data Release Tag
^^^^^^^^^^^^^^^^

.. code-block:: yaml

   release: npipe6v20

String identifier included in output filenames. Use the official Planck
release name for consistency with published data.

Computational Settings
^^^^^^^^^^^^^^^^^^^^^^

``test``
   **Test mode** (boolean, default: ``false``)
   
   * ``true``: Sparse pixel sampling (every 64th pixel) for fast debugging
   * ``false``: Full pixel sampling for production runs

``conserve_memory``
   **Memory conservation mode** (boolean, default: ``true``)
   
   * ``true``: Process hit maps in smaller chunks (reduces RAM usage)
   * ``false``: Keep full hit maps in memory (faster but uses more RAM)
   
   Recommended ``true`` for large detector sets or limited RAM systems.

``overwrite``
   **Overwrite existing files** (boolean, default: ``false``)
   
   * ``true``: Recompute and overwrite existing NPZ and FITS files
   * ``false``: Skip detector pairs that already have output files

``force_det``
   **Force detector label** (string or ``null``, default: ``null``)
   
   If set, adds ``_FD<label>`` to filenames. Used for debugging or special runs.

Output File Options
^^^^^^^^^^^^^^^^^^^

Control which FITS output files are generated:

.. code-block:: yaml

   blfile: true       # Scalar B_l(T) window function
   blTEBfile: true    # Separate B_l for T, E, B
   wlfile: true       # Full W_l coupling matrices
   do_plot: false     # Generate PNG plots of beam windows

Set any to ``false`` to skip that output type.

Detector Pairs
^^^^^^^^^^^^^^

List of ``[detset1, detset2]`` pairs to process:

.. code-block:: yaml

   detpairs:
     # Auto-spectra
     - ['143A', '143A']
     - ['143B', '143B']
     
     # Cross-spectra
     - ['143A', '143B']
     
     # Single detectors
     - ['143-1a', '143-1a']
     - ['143-1a', '143-1b']
     
     # Frequency cross-spectra
     - ['143GHz', '217GHz']

Valid detector set names:

* Full frequency: ``30GHz``, ``44GHz``, ``70GHz``, ``100GHz``, ``143GHz``, 
  ``217GHz``, ``353GHz``, ``545GHz``, ``857GHz``
* Subset A/B: ``100A``, ``100B``, ``143A``, ``143B``, etc.
* Single detectors: ``LFI27M``, ``143-1a``, etc.
* Single horns: ``100-1``, ``143-2`` (combines ``a`` and ``b``)

.. tip::
   Use the helper script ``scripts/make_planck_detpairs.py`` to generate
   the full list of detector pairs used in the original NPIPE processing.

Example Configurations
----------------------

Minimal Configuration
^^^^^^^^^^^^^^^^^^^^^

For a quick test with a single detector pair:

.. code-block:: yaml

   rimo_hfi: /data/RIMO_HFI.fits
   blmfile: /data/beams/blm_{}.fits
   outdir: ./test_output
   
   smax: 6
   release: test
   test: true
   overwrite: true
   
   detpairs:
     - ['143-1a', '143-1a']

Full Production Run
^^^^^^^^^^^^^^^^^^^^

Complete configuration for Planck 143 GHz analysis:

.. code-block:: yaml

   # RIMO
   rimo_lfi: /data/planck/RIMO_LFI_npipe5_symmetrized.fits
   rimo_hfi: /data/planck/RIMO_HFI_npipe5v16_symmetrized.fits
   
   # Beams
   blmfile: /data/planck/beams/blm_{}.fits
   
   # I/O
   outdir: /scratch/beam_windows/npipe6v20
   indir: /scratch/beam_windows/npipe6v20
   
   # Physics
   smax: 6
   spin_ref: Pxx
   blm_ref: Dxx
   angle_shift: 0.0
   
   # Models
   release: npipe6v20
   rhobeam: IMO
   rhohit: IMO
   
   # Computational
   test: false
   conserve_memory: true
   overwrite: false
   
   # Outputs
   blfile: true
   blTEBfile: true
   wlfile: true
   do_plot: true
   
   # Detector pairs
   detpairs:
     # Auto-spectra
     - ['143GHz', '143GHz']
     - ['143A', '143A']
     - ['143B', '143B']
     # Cross-spectra
     - ['143A', '143B']
     # All individual detectors
     - ['143-1a', '143-1a']
     - ['143-1b', '143-1b']
     - ['143-2a', '143-2a']
     - ['143-2b', '143-2b']
     - ['143-3a', '143-3a']
     - ['143-3b', '143-3b']
     - ['143-4a', '143-4a']
     - ['143-4b', '143-4b']
     - ['143-5', '143-5']
     - ['143-6', '143-6']
     - ['143-7', '143-7']

Multi-Frequency Cross-Spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For component separation analyses:

.. code-block:: yaml

   rimo_hfi: /data/RIMO_HFI.fits
   blmfile: /data/beams/blm_{}.fits
   outdir: ./freq_cross_beams
   
   smax: 6
   release: npipe6v20
   
   detpairs:
     # All frequency auto-spectra
     - ['100GHz', '100GHz']
     - ['143GHz', '143GHz']
     - ['217GHz', '217GHz']
     - ['353GHz', '353GHz']
     
     # All cross-spectra
     - ['100GHz', '143GHz']
     - ['100GHz', '217GHz']
     - ['143GHz', '217GHz']
     - ['143GHz', '353GHz']
     - ['217GHz', '353GHz']

Running from YAML
-----------------

Command line:

.. code-block:: bash

   python -m qp_planck.qp_pipeline config.yaml

From Python:

.. code-block:: python

   from qp_planck import run_qp_from_yaml
   
   run_qp_from_yaml("config.yaml")

With MPI:

.. code-block:: bash

   mpirun -n 16 python -m qp_planck.qp_pipeline config.yaml

Validation
----------

The pipeline will validate your configuration and report errors for:

* Missing RIMO files
* Invalid detector names
* Incompatible parameter combinations
* Missing input data files

Check the console output for warnings and error messages.

Best Practices
--------------

1. **Start small**: Test with ``test: true`` and a single detector pair first
2. **Organize outputs**: Use descriptive ``release`` tags and separate ``outdir`` paths
3. **Version control**: Keep YAML configs in git alongside analysis scripts
4. **Document runs**: Add comments to YAML files explaining non-standard choices
5. **Reproducibility**: Record exact paths and file versions for published results

See Also
--------

* :doc:`quickstart` — Basic usage examples
* :doc:`examples` — Advanced configuration scenarios
* :doc:`api/pipeline` — API reference for configuration options
