Configuration Guide
===================

This page details the configuration options for **qp-planck** when using
the Python API directly (without YAML).

Function Signature
------------------

.. code-block:: python

   from qp_planck import run_qp_pipeline
   
   run_qp_pipeline(
       detpairs,
       rimo=None,
       rimo_lfi=None,
       rimo_hfi=None,
       blmfile="blm_{}.fits",
       outdir="../quickpol_output",
       indir=None,
       smax=6,
       spin_ref="Pxx",
       blm_ref="Dxx",
       angle_shift=0.0,
       force_det=None,
       release="npipe6v20",
       rhobeam="IMO",
       rhohit="IMO",
       test=False,
       planet="",
       conserve_memory=True,
       overwrite=False,
       blfile=True,
       blTEBfile=True,
       wlfile=True,
       do_plot=False,
       config=None
   )

Required Parameters
-------------------

``detpairs``
   **Detector pair list** (list of tuples)
   
   List of ``(detset1, detset2)`` tuples specifying which detector combinations
   to process.
   
   Example::
   
      detpairs = [
          ("143A", "143A"),      # Auto-spectrum
          ("143A", "143B"),      # Cross-spectrum
          ("143-1a", "143-1a"),  # Single detector
      ]
   
   See :doc:`yaml_config` for valid detector set names.

RIMO Configuration
------------------

Provide at least one of:

``rimo``
   **Single merged RIMO file** (string, default: None)
   
   Path to a FITS file containing the full Planck RIMO (LFI + HFI merged).

``rimo_lfi``
   **LFI RIMO file** (string, default: None)
   
   Path to the LFI Reduced Instrument Model FITS file.

``rimo_hfi``
   **HFI RIMO file** (string, default: None)
   
   Path to the HFI Reduced Instrument Model FITS file.

Example::

   # Option 1: Separate files
   run_qp_pipeline(
       detpairs=pairs,
       rimo_lfi="/data/RIMO_LFI_npipe5.fits",
       rimo_hfi="/data/RIMO_HFI_npipe5v16.fits",
       ...
   )
   
   # Option 2: Merged file
   run_qp_pipeline(
       detpairs=pairs,
       rimo="/data/RIMO_Planck_merged.fits",
       ...
   )

Beam and I/O Paths
------------------

``blmfile``
   **Beam multipole file template** (string, default: ``"blm_{}.fits"``)
   
   Path template for beam multipole FITS files. The ``{}`` is replaced by
   the detector name.
   
   Example::
   
      blmfile="/data/planck/beams/blm_{}.fits"
      # Expands to: blm_143-1a.fits, blm_143-1b.fits, etc.
   
   Set to ``""`` to use Gaussian beams from RIMO FWHM instead.

``outdir``
   **Output directory** (string, default: ``"../quickpol_output"``)
   
   Directory where NPZ beam matrices and FITS window functions will be written.
   Created automatically if it doesn't exist.

``indir``
   **Input directory for NPZ files** (string or None, default: None)
   
   Directory to read NPZ beam matrices from when running the mat2fits step.
   If ``None``, defaults to ``outdir``.

Physics Parameters
------------------

``smax``
   **Maximum spin moment** (int, default: 6)
   
   Maximum spin :math:`s` in the spin-weighted harmonic expansion. Values
   0–6 are stored in hit maps. Higher values cost more memory but are
   typically not needed for Planck-level systematics.

``spin_ref``
   **Spin moment reference frame** (``'Pxx'`` or ``'Dxx'``, default: ``'Pxx'``)
   
   Reference frame for the spin moments in the hit maps:
   
   * ``'Pxx'``: Detector polarization-sensitive direction
   * ``'Dxx'``: Detector :math:`uv` frame
   
   Must match the convention used when generating the polmoments files.

``blm_ref``
   **Beam reference frame** (``'Pxx'`` or ``'Dxx'``, default: ``'Dxx'``)
   
   Reference frame in which the beam multipoles :math:`b_{\ell m}` are defined.

``angle_shift``
   **Detector angle shift** (float, degrees, default: 0.0)
   
   Additional rotation angle applied to all detectors. Used for testing
   systematic uncertainties in detector orientation.

``planet``
   **Planet deconvolution** (string, default: ``""``)
   
   Name of planet to deconvolve from beam (e.g., ``"Saturn"``). Empty string
   means no planetary deconvolution. Experimental feature.

Polarization Models
-------------------

``rhobeam``
   **Beam polarization model** (``'IMO'``, ``'Ideal'``, or float, default: ``'IMO'``)
   
   Model for the beam cross-polarization response :math:`\rho`:
   
   * ``'IMO'``: Use :math:`\rho = (1-\epsilon)/(1+\epsilon)` from RIMO
   * ``'Ideal'``: :math:`\rho = 1` for PSBs, :math:`\rho = 0` for SWBs
   * Float value: Fixed :math:`\rho` for all detectors (debugging)

``rhohit``
   **Hit matrix polarization model** (``'IMO'``, ``'Ideal'``, or float, default: ``'IMO'``)
   
   Same as ``rhobeam`` but applied to the hit matrix calculation. Usually
   set to the same value as ``rhobeam``.

Data Release
------------

``release``
   **Release tag** (string, default: ``"npipe6v20"``)
   
   String identifier for the data release. Included in output filenames.
   Use official Planck release names for consistency.

``force_det``
   **Force detector label** (string or None, default: None)
   
   If set, adds ``_FD<label>`` to output filenames. For debugging or
   special processing runs.

Computational Options
---------------------

``test``
   **Test mode** (bool, default: False)
   
   * ``True``: Use sparse pixel sampling (every 64th pixel) for fast testing
   * ``False``: Full pixel sampling for production runs

``conserve_memory``
   **Memory conservation** (bool, default: True)
   
   * ``True``: Process hit maps in chunks to reduce RAM usage
   * ``False``: Keep full hit maps in memory (faster but uses more RAM)

``overwrite``
   **Overwrite existing files** (bool, default: False)
   
   * ``True``: Recompute and overwrite existing NPZ and FITS files
   * ``False``: Skip detector pairs with existing outputs

Output File Options
-------------------

``blfile``
   **Generate scalar B_l FITS** (bool, default: True)
   
   Create temperature beam window function ``Bl_<release>_<det1>x<det2>.fits``.

``blTEBfile``
   **Generate TEB B_l FITS** (bool, default: True)
   
   Create separate T, E, B beam windows ``Bl_TEB_<release>_<det1>x<det2>.fits``.

``wlfile``
   **Generate W_l matrix FITS** (bool, default: True)
   
   Create full coupling matrix ``Wl_<release>_<det1>x<det2>.fits``.

``do_plot``
   **Generate plots** (bool, default: False)
   
   If ``True``, also create PNG plots of the beam window functions.
   Requires matplotlib.

YAML Override
-------------

``config``
   **YAML configuration file or dict** (string, Path, or dict, default: None)
   
   If provided, values from the YAML file (or dict) are used as defaults,
   but explicit keyword arguments take precedence.
   
   Example::
   
      run_qp_pipeline(
          detpairs=[],  # will be read from YAML
          config="beam_config.yaml",
          overwrite=True  # overrides YAML value
      )

Complete Example
----------------

Full configuration for a production run:

.. code-block:: python

   from qp_planck import run_qp_pipeline
   
   # Define detector pairs
   pairs = [
       ("143GHz", "143GHz"),
       ("143A", "143A"),
       ("143B", "143B"),
       ("143A", "143B"),
   ]
   
   # Run pipeline
   run_qp_pipeline(
       detpairs=pairs,
       
       # RIMO
       rimo_lfi="/data/planck/RIMO_LFI_npipe5_symmetrized.fits",
       rimo_hfi="/data/planck/RIMO_HFI_npipe5v16_symmetrized.fits",
       
       # Beams
       blmfile="/data/planck/beams/blm_{}.fits",
       
       # I/O
       outdir="/scratch/beam_windows/npipe6v20",
       indir="/scratch/beam_windows/npipe6v20",
       
       # Physics
       smax=6,
       spin_ref="Pxx",
       blm_ref="Dxx",
       angle_shift=0.0,
       planet="",
       
       # Models
       release="npipe6v20",
       rhobeam="IMO",
       rhohit="IMO",
       
       # Computational
       test=False,
       conserve_memory=True,
       overwrite=False,
       force_det=None,
       
       # Outputs
       blfile=True,
       blTEBfile=True,
       wlfile=True,
       do_plot=True,
   )

Hybrid Configuration
--------------------

You can combine YAML and Python arguments:

.. code-block:: python

   from qp_planck import run_qp_pipeline
   
   # Load base config from YAML, override specific parameters
   run_qp_pipeline(
       detpairs=[],  # read from YAML
       config="base_config.yaml",
       
       # Override for this run
       test=True,
       overwrite=True,
       outdir="./debug_output"
   )

This is useful for:

* Testing with different settings
* Running on different systems (different paths)
* Quick parameter scans

See Also
--------

* :doc:`yaml_config` — YAML configuration reference
* :doc:`api/pipeline` — API documentation for ``run_qp_pipeline``
* :doc:`examples` — Configuration examples for common use cases
