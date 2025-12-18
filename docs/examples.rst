Examples
========

This page provides complete working examples for common use cases.

Example 1: Single Detector Auto-Spectrum
-----------------------------------------

Compute the beam window function for a single Planck detector:

.. code-block:: python

   from qp_planck import run_qp_pipeline
   
   # Paths to data files
   RIMO_HFI = "/data/planck/RIMO_HFI_npipe5v16_symmetrized.fits"
   BEAMS = "/data/planck/beams/blm_{}.fits"
   OUTPUT = "./detector_beams"
   
   # Process detector 143-1a
   run_qp_pipeline(
       detpairs=[("143-1a", "143-1a")],
       rimo_hfi=RIMO_HFI,
       blmfile=BEAMS,
       outdir=OUTPUT,
       smax=6,
       spin_ref="Pxx",
       blm_ref="Dxx",
       release="npipe6v20",
       test=False,
       overwrite=True
   )
   
   # Results:
   # - detector_beams/beam_matrix_143-1ax143-1a_*.npz
   # - detector_beams/Bl_npipe6v20_143-1ax143-1a.fits
   # - detector_beams/Bl_TEB_npipe6v20_143-1ax143-1a.fits
   # - detector_beams/Wl_npipe6v20_143-1ax143-1a.fits

Example 2: Frequency Channel Cross-Spectra
-------------------------------------------

Compute all auto and cross-spectra for two frequency channels:

.. code-block:: python

   from qp_planck import run_qp_pipeline
   
   # All 143x217 combinations
   pairs = [
       ("143GHz", "143GHz"),  # 143 auto
       ("217GHz", "217GHz"),  # 217 auto
       ("143GHz", "217GHz"),  # 143x217 cross
   ]
   
   run_qp_pipeline(
       detpairs=pairs,
       rimo_hfi="/data/RIMO_HFI.fits",
       blmfile="/data/beams/blm_{}.fits",
       outdir="./freq_cross_beams",
       release="npipe6v20_143x217",
       conserve_memory=True,
       overwrite=False
   )

Example 3: Detector Subset Comparison
--------------------------------------

Compare subset A and B within a frequency:

.. code-block:: python

   from qp_planck import run_qp_pipeline
   
   # 143 GHz subsets
   pairs = [
       ("143A", "143A"),  # Subset A auto
       ("143B", "143B"),  # Subset B auto
       ("143A", "143B"),  # A x B cross
   ]
   
   run_qp_pipeline(
       detpairs=pairs,
       rimo_hfi="/data/RIMO_HFI.fits",
       blmfile="/data/beams/blm_{}.fits",
       outdir="./subset_comparison",
       release="npipe6v20_subsets"
   )
   
   # Use the cross-spectrum to test for systematic differences

Example 4: YAML-Based Production Run
-------------------------------------

Large-scale processing with YAML configuration:

**config_production.yaml:**

.. code-block:: yaml

   # Full production configuration
   rimo_lfi: /data/planck/RIMO_LFI_npipe5_symmetrized.fits
   rimo_hfi: /data/planck/RIMO_HFI_npipe5v16_symmetrized.fits
   blmfile: /data/planck/beams/blm_{}.fits
   
   outdir: /scratch/beam_windows/npipe6v20
   
   smax: 6
   spin_ref: Pxx
   blm_ref: Dxx
   release: npipe6v20
   rhobeam: IMO
   rhohit: IMO
   
   test: false
   conserve_memory: true
   overwrite: false
   
   blfile: true
   blTEBfile: true
   wlfile: true
   do_plot: true
   
   detpairs:
     # Full frequency auto-spectra
     - ['100GHz', '100GHz']
     - ['143GHz', '143GHz']
     - ['217GHz', '217GHz']
     - ['353GHz', '353GHz']
     
     # Subset auto-spectra for CMB channels
     - ['100A', '100A']
     - ['100B', '100B']
     - ['143A', '143A']
     - ['143B', '143B']
     - ['217A', '217A']
     - ['217B', '217B']
     
     # Cross-spectra
     - ['100A', '100B']
     - ['143A', '143B']
     - ['217A', '217B']

**Run with Python:**

.. code-block:: python

   from qp_planck import run_qp_from_yaml
   
   run_qp_from_yaml("config_production.yaml")

**Run with MPI:**

.. code-block:: bash

   mpirun -n 32 python -m qp_planck.qp_pipeline config_production.yaml

Example 5: Gaussian Beam Test
------------------------------

Use Gaussian beams instead of measured beam multipoles:

.. code-block:: python

   from qp_planck import run_qp_pipeline
   
   # Empty blmfile = use Gaussian beams from RIMO FWHM
   run_qp_pipeline(
       detpairs=[("143A", "143A")],
       rimo_hfi="/data/RIMO_HFI.fits",
       blmfile="",  # Empty = Gaussian
       outdir="./gaussian_beams",
       release="gaussian_test",
       overwrite=True
   )

Example 6: Systematic Tests
----------------------------

Test impact of detector angle uncertainty:

.. code-block:: python

   from qp_planck import run_qp_pipeline
   import numpy as np
   
   angles = [0.0, 0.1, 0.2, 0.5, 1.0]  # degrees
   
   for angle in angles:
       run_qp_pipeline(
           detpairs=[("143A", "143A")],
           rimo_hfi="/data/RIMO_HFI.fits",
           blmfile="/data/beams/blm_{}.fits",
           outdir=f"./angle_systematic",
           angle_shift=angle,
           release=f"angle_{angle:.1f}deg",
           overwrite=True
       )
   
   # Compare outputs to assess sensitivity to angle errors

Example 7: Analyzing Outputs
-----------------------------

Complete analysis workflow:

.. code-block:: python

   from astropy.io import fits
   import matplotlib.pyplot as plt
   import numpy as np
   
   # Read beam windows
   def read_beam_windows(filename):
       with fits.open(filename) as hdul:
           bl_T = hdul[1].data['T']
           bl_E = hdul[1].data['E']
           bl_B = hdul[1].data['B']
       return bl_T, bl_E, bl_B
   
   # Compare detector pairs
   bl_1a_T, bl_1a_E, bl_1a_B = read_beam_windows(
       "Bl_TEB_npipe6v20_143-1ax143-1a.fits"
   )
   bl_1b_T, bl_1b_E, bl_1b_B = read_beam_windows(
       "Bl_TEB_npipe6v20_143-1bx143-1b.fits"
   )
   
   ell = np.arange(len(bl_1a_T))
   
   # Plot comparison
   fig, axes = plt.subplots(1, 3, figsize=(15, 5))
   
   for ax, bl_a, bl_b, label in zip(
       axes,
       [bl_1a_T, bl_1a_E, bl_1a_B],
       [bl_1b_T, bl_1b_E, bl_1b_B],
       ['T', 'E', 'B']
   ):
       ax.semilogy(ell, bl_a, label='143-1a')
       ax.semilogy(ell, bl_b, label='143-1b')
       ax.set_xlabel(r'$\ell$')
       ax.set_ylabel(r'$B_\ell$')
       ax.set_title(f'{label} Beam Windows')
       ax.legend()
       ax.grid(alpha=0.3)
       ax.set_xlim(0, 2000)
       ax.set_ylim(1e-8, 1)
   
   plt.tight_layout()
   plt.savefig('beam_comparison_143.png', dpi=150)
   plt.show()
   
   # Compute fractional difference
   diff_T = (bl_1a_T - bl_1b_T) / bl_1a_T
   
   plt.figure()
   plt.plot(ell, diff_T * 100)
   plt.xlabel(r'$\ell$')
   plt.ylabel('Fractional Difference (%)')
   plt.title('143-1a vs 143-1b Beam Mismatch')
   plt.axhline(0, color='k', linestyle='--')
   plt.grid(alpha=0.3)
   plt.xlim(0, 2000)
   plt.savefig('beam_mismatch_143.png', dpi=150)
   plt.show()

Example 8: Applying Beam Windows
---------------------------------

Use beam windows in power spectrum analysis:

.. code-block:: python

   import healpy as hp
   import numpy as np
   from astropy.io import fits
   
   # Load CMB map (temperature)
   cmb_map = hp.read_map("cmb_143GHz_map.fits")
   
   # Compute angular power spectrum
   cl_obs = hp.anafast(cmb_map)
   
   # Load beam window
   with fits.open("Bl_npipe6v20_143GHzx143GHz.fits") as hdul:
       bl = hdul[1].data['TEMPERATURE']
   
   # Deconvolve beam
   cl_deconv = cl_obs[:len(bl)] / bl**2
   
   # Handle zeros in bl (high ℓ)
   valid = bl > 1e-10
   cl_deconv[~valid] = 0
   
   # Plot
   ell = np.arange(len(cl_obs))
   
   plt.figure(figsize=(10, 6))
   plt.loglog(ell[2:], ell[2:] * (ell[2:] + 1) * cl_obs[2:] / (2*np.pi),
              label='Observed')
   plt.loglog(ell[2:], ell[2:] * (ell[2:] + 1) * cl_deconv[2:] / (2*np.pi),
              label='Deconvolved')
   plt.xlabel(r'$\ell$')
   plt.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$ [$\mu$K$^2$]')
   plt.title('143 GHz Power Spectrum')
   plt.legend()
   plt.grid(alpha=0.3)
   plt.xlim(2, 2000)
   plt.savefig('power_spectrum_143.png', dpi=150)
   plt.show()

Example 9: Memory-Efficient Processing
---------------------------------------

Process large detector sets on limited-RAM systems:

.. code-block:: python

   from qp_planck import run_qp_pipeline, list_planck
   
   # Get all individual detectors for 143 GHz
   dets_143 = list_planck("143GHz")
   
   # Create auto-spectrum pairs
   pairs = [(det, det) for det in dets_143]
   
   print(f"Processing {len(pairs)} detector auto-spectra...")
   
   # Process with memory conservation
   run_qp_pipeline(
       detpairs=pairs,
       rimo_hfi="/data/RIMO_HFI.fits",
       blmfile="/data/beams/blm_{}.fits",
       outdir="./individual_detectors",
       conserve_memory=True,  # Use less RAM
       test=False,
       release="npipe6v20_individuals"
   )

Example 10: HPC Batch Script
-----------------------------

SLURM script for large-scale HPC processing:

**run_beams.sh:**

.. code-block:: bash

   #!/bin/bash
   #SBATCH --job-name=qp_beams
   #SBATCH --nodes=4
   #SBATCH --ntasks=128
   #SBATCH --time=24:00:00
   #SBATCH --mem-per-cpu=4G
   #SBATCH --output=qp_beams_%j.log
   
   module load python/3.11
   module load openmpi/4.1
   
   # Activate environment
   source /path/to/venv/bin/activate
   
   # Run pipeline with MPI
   mpirun -n $SLURM_NTASKS python -m qp_planck.qp_pipeline \
       config_full_planck.yaml
   
   echo "Job finished: $(date)"

Submit with::

   sbatch run_beams.sh

See Also
--------

* :doc:`quickstart` — Basic introduction
* :doc:`yaml_config` — Configuration reference
* :doc:`api/pipeline` — API documentation
