API Reference: hmap2mat
========================

.. currentmodule:: qp_planck.qp_hmap2mat

Low-level QuickPol driver for building beam matrices from hitmaps and beam
multipoles.

Main Function
-------------

.. autofunction:: hmap2mat

Core Beam Matrix Functions
---------------------------

.. autofunction:: fill_beam_dict

.. autofunction:: bmat

.. autofunction:: make_hit_matrix

.. autofunction:: product_pre2

Hit Matrix Operations
---------------------

.. autofunction:: make_hit_vectors

.. autofunction:: invert_hit

.. autofunction:: invert_hit_sub

.. autofunction:: apply_hit_masks

Beam Reading Functions
----------------------

.. autofunction:: get_blm_det

.. autofunction:: get_gblm_det

.. autofunction:: get_wsmap_det

Matrix Operations
-----------------

.. autofunction:: dots

.. autofunction:: adjoint

.. autofunction:: interpolate_matrix

.. autofunction:: deconv_planet

Helper Functions
----------------

.. autofunction:: parse_detname

.. autofunction:: get_angles

.. autofunction:: my_create_bololist_w8

.. autofunction:: proc_angle

.. autofunction:: count_pix

.. autofunction:: get_all_masks

.. autofunction:: print_matrix

Utility Functions
-----------------

.. autofunction:: detset2nside

Constants
---------

This module defines several color codes for terminal output:

.. code-block:: python

   NO_COLOR = "\\x1b[0m"
   RED_COLOR = "\\x1b[31;01m"
   GREEN_COLOR = "\\x1b[32;11m"
   YELLOW_COLOR = "\\x1b[33;11m"
   BLUE_COLOR = "\\x1b[34;11m"
   MAGENTA_COLOR = "\\x1b[35;11m"
   CYAN_COLOR = "\\x1b[36;11m"

Examples
--------

Direct hmap2mat Usage
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck import load_RIMO, hmap2mat
   
   # Load RIMO
   RIMO_HFI = load_RIMO("RIMO_HFI_npipe5v16.fits")
   RIMO_LFI = load_RIMO("RIMO_LFI_npipe5.fits")
   RIMO = {**RIMO_LFI, **RIMO_HFI}
   
   # Compute beam matrix for a detector pair
   hmap2mat(
       RIMO=RIMO,
       mytype=("143A", "143B"),  # detector pair
       blmfile="/data/beams/blm_{}.fits",
       outdir="./beam_matrices",
       smax=6,
       spin_ref="Pxx",
       blm_ref="Dxx",
       angle_shift=0.0,
       release="npipe6v20",
       rhobeam="IMO",
       rhohit="IMO",
       test=False,
       conserve_memory=True,
       overwrite=False
   )

Working with Beam Multipoles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck.qp_hmap2mat import get_blm_det, get_gblm_det
   
   # Read measured beam multipoles
   blm = get_blm_det(
       blmfile="/data/beams/blm_{}.fits",
       det="143-1a",
       lmax=4000,
       mmax=10
   )
   print(f"Beam shape: {blm.shape}")  # (lmax+1, mmax+1, ndb)
   
   # Generate Gaussian beam from RIMO FWHM
   from qp_planck import load_RIMO
   RIMO = load_RIMO("RIMO_HFI.fits")
   
   blm_gauss = get_gblm_det(
       RIMO=RIMO,
       det="143-1a",
       lmax=4000,
       mmax=0  # Gaussian is m=0 only
   )

Constructing Beam Matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck import load_RIMO
   from qp_planck.qp_hmap2mat import fill_beam_dict, bmat
   
   # Load RIMO
   RIMO = load_RIMO("RIMO_HFI.fits")
   
   # Build beam dictionaries for all detectors in a set
   bdict = fill_beam_dict(
       RIMO=RIMO,
       blmfile="/data/beams/blm_{}.fits",
       lmax=4000,
       mmax=10,
       detset="143A",
       blm_ref="Dxx",
       angle_sdeg=0.0
   )
   
   # Construct beam matrix for specific ℓ and spin s
   ell = 100
   spin = 0
   bm = bmat(bdict[0], ell, spin, rhobeam="IMO", verbose=True)
   print(f"Beam matrix at ℓ={ell}, s={spin}:")
   print(bm)

Hit Matrix Computation
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck import load_RIMO
   from qp_planck.qp_hmap2mat import make_hit_matrix
   
   RIMO = load_RIMO("RIMO_HFI.fits")
   
   # Compute hit matrix for detector pair
   hitmat, varmat, v2mat, nbad1, nbad2, skip = make_hit_matrix(
       RIMO=RIMO,
       nside=2048,
       hgrp1="/data/polmoments",
       detset1="143A",
       hgrp2="/data/polmoments",
       detset2="143A",
       smax=6,
       spin_ref="Pxx",
       test=False,
       thr=3e-3,
       angle_shift=0.0,
       savefile="",
       release="npipe6v20",
       rhohit="IMO",
       masks=[None, None],
       conserve_memory=True
   )
   
   print(f"Hit matrix shape: {hitmat.shape}")
   print(f"Bad pixels: {nbad1}, {nbad2}")

Advanced: Custom Masking
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from qp_planck.qp_hmap2mat import get_all_masks
   import healpy as hp
   
   # Load custom masks
   mask_T = hp.read_map("mask_temperature.fits")
   mask_P = hp.read_map("mask_polarization.fits")
   
   # Apply masks in hmap2mat
   # (Currently get_all_masks is a stub returning no masks,
   #  but you can modify it for custom masking)
   masks, mask_names, w_cutsky, mask_means = get_all_masks(
       release="custom",
       dets=["143A", "143A"],
       lmax=4000
   )
