API Reference: mat2fits
========================

.. currentmodule:: qp_planck.qp_mat2fits

Converter from QuickPol NPZ beam matrices to FITS window function files.

Main Function
-------------

.. autofunction:: mat2fits

FITS Writing Functions
----------------------

.. autofunction:: my_mwrfits

.. autofunction:: clobber

Beam Fitting
------------

.. autofunction:: fit_gauss

Helper Functions
----------------

.. autofunction:: detset2lmax

.. autofunction:: detset2pol

Examples
--------

Converting NPZ to FITS
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck.qp_mat2fits import mat2fits
   
   # Convert beam matrix NPZ to FITS window functions
   mat2fits(
       indir="./beam_matrices",
       outdir="./beam_windows",
       dets=("143A", "143A"),
       smax=6,
       release="npipe6v20",
       full=True,
       blfile=True,      # Create Bl_*.fits
       blTEBfile=True,   # Create Bl_TEB_*.fits
       wlfile=True,      # Create Wl_*.fits
       overwrite=True,
       do_plot=True      # Also make PNG plots
   )

This produces:

* ``Bl_npipe6v20_143Ax143A.fits`` — Scalar :math:`B_\ell(T)`
* ``Bl_TEB_npipe6v20_143Ax143A.fits`` — :math:`B_\ell` for T, E, B
* ``Wl_npipe6v20_143Ax143A.fits`` — Full :math:`W_\ell` matrices
* ``Bl_npipe6v20_143Ax143A.png``, etc. — Plots (if ``do_plot=True``)

Reading FITS Outputs
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from astropy.io import fits
   import matplotlib.pyplot as plt
   import numpy as np
   
   # Read scalar B_l
   with fits.open("Bl_npipe6v20_143Ax143A.fits") as hdul:
       bl_T = hdul[1].data['TEMPERATURE']
   
   # Read TEB B_l
   with fits.open("Bl_TEB_npipe6v20_143Ax143A.fits") as hdul:
       bl_T = hdul[1].data['T']
       bl_E = hdul[1].data['E']
       bl_B = hdul[1].data['B']
   
   # Read W_l matrices
   with fits.open("Wl_npipe6v20_143Ax143A.fits") as hdul:
       # Extension 1: TT coupling
       wl_TT_to_TT = hdul[1].data['TT_2_TT']
       wl_TT_to_TE = hdul[1].data['TT_2_TE']
       
       # Extension 2: EE coupling
       wl_EE_to_EE = hdul[2].data['EE_2_EE']
       wl_EE_to_BB = hdul[2].data['EE_2_BB']
       
       # Extension 3: BB coupling
       wl_BB_to_BB = hdul[3].data['BB_2_BB']
       
       # Extension 4: TE coupling
       wl_TE_to_TE = hdul[4].data['TE_2_TE']
   
   # Plot beam windows
   ell = np.arange(len(bl_T))
   plt.figure(figsize=(10, 6))
   plt.semilogy(ell, bl_T, label='T')
   plt.semilogy(ell, bl_E, label='E')
   plt.semilogy(ell, bl_B, label='B')
   plt.xlabel(r'$\ell$')
   plt.ylabel(r'$B_\ell$')
   plt.legend()
   plt.xlim(0, 3000)
   plt.ylim(1e-8, 1)
   plt.grid(alpha=0.3)
   plt.title('143A Beam Windows')
   plt.show()

Fitting Gaussian Beams
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck.qp_mat2fits import fit_gauss
   from astropy.io import fits
   import numpy as np
   
   # Read measured beam window
   with fits.open("Bl_npipe6v20_143-1ax143-1a.fits") as hdul:
       bl_measured = hdul[1].data['TEMPERATURE']
   
   # Fit Gaussian
   bl_gauss, sigma_rad = fit_gauss(bl_measured)
   
   # Convert to FWHM in arcminutes
   fwhm_arcmin = np.degrees(sigma_rad) * 60 * np.sqrt(8 * np.log(2))
   print(f"Best-fit Gaussian FWHM: {fwhm_arcmin:.2f} arcmin")
   
   # Plot comparison
   ell = np.arange(len(bl_measured))
   plt.semilogy(ell, bl_measured, label='Measured')
   plt.semilogy(ell, bl_gauss, '--', label=f'Gaussian {fwhm_arcmin:.2f}\'')
   plt.xlabel(r'$\ell$')
   plt.ylabel(r'$B_\ell$')
   plt.legend()
   plt.show()

Batch Processing
^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck.qp_mat2fits import mat2fits
   
   # Process all 143 GHz detector pairs
   pairs = [
       ("143-1a", "143-1a"),
       ("143-1b", "143-1b"),
       ("143-2a", "143-2a"),
       ("143-2b", "143-2b"),
       # ... etc
   ]
   
   for det1, det2 in pairs:
       print(f"Processing {det1} x {det2}...")
       mat2fits(
           indir="./beam_matrices",
           outdir="./beam_windows",
           dets=(det1, det2),
           smax=6,
           release="npipe6v20",
           full=True,
           overwrite=False  # Skip if already exists
       )

Custom FITS Writing
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck.qp_mat2fits import my_mwrfits
   import numpy as np
   
   # Create custom beam window FITS file
   lmax = 3000
   ell = np.arange(lmax + 1)
   
   # Synthetic Gaussian beam
   fwhm_rad = np.radians(7.0 / 60)  # 7 arcmin
   bl = np.exp(-0.5 * ell * (ell + 1) * fwhm_rad**2)
   
   # Write to FITS
   my_mwrfits(
       filename="custom_beam.fits",
       data=[[bl]],
       colnames=[["TEMPERATURE"]],
       bintable=False,
       ftype="B",
       extnames=["WINDOW FUNCTION"],
       origin=["Custom synthetic beam"],
       dets=["CUSTOM-1"]
   )

Checking Outputs
^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck.qp_mat2fits import clobber
   
   filename = "Bl_npipe6v20_143Ax143A.fits"
   
   # Check if file exists and whether to write
   should_write = clobber(filename, overwrite=False)
   
   if should_write:
       print(f"Will write {filename}")
   else:
       print(f"{filename} exists, skipping")
