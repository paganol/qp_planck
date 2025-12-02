References
==========

Scientific Papers
-----------------

QuickPol Method
^^^^^^^^^^^^^^^

.. [Hivon2017] Hivon, E., Mottet, S., & Ponthieu, N. 2017,
   "QuickPol: Fast calculation of effective beam matrices for CMB polarization",
   *Astronomy & Astrophysics*, **598**, A25
   
   DOI: `10.1051/0004-6361/201629626 <https://doi.org/10.1051/0004-6361/201629626>`_
   
   ArXiv: `1606.09313 <https://arxiv.org/abs/1606.09313>`_

Planck Mission and Data Releases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. [Planck2020_LVII] Planck Collaboration Int. LVII. 2020,
   "Planck intermediate results. LVII. Joint Planck LFI and HFI data processing",
   *Astronomy & Astrophysics*, **643**, A42
   
   DOI: `10.1051/0004-6361/202038073 <https://doi.org/10.1051/0004-6361/202038073>`_

.. [Planck2016_IV] Planck Collaboration IV. 2016,
   "Planck 2015 results. IV. Low Frequency Instrument beams and window functions",
   *Astronomy & Astrophysics*, **594**, A4
   
   DOI: `10.1051/0004-6361/201525809 <https://doi.org/10.1051/0004-6361/201525809>`_

.. [Planck2016_VII] Planck Collaboration VII. 2016,
   "Planck 2015 results. VII. High Frequency Instrument data processing: Time-ordered information and beams",
   *Astronomy & Astrophysics*, **594**, A7
   
   DOI: `10.1051/0004-6361/201525844 <https://doi.org/10.1051/0004-6361/201525844>`_

.. [Planck2016_I] Planck Collaboration I. 2016,
   "Planck 2015 results. I. Overview of products and scientific results",
   *Astronomy & Astrophysics*, **594**, A1
   
   DOI: `10.1051/0004-6361/201527101 <https://doi.org/10.1051/0004-6361/201527101>`_

HEALPix and Spherical Harmonics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. [Gorski2005] Górski, K. M., Hivon, E., Banday, A. J., et al. 2005,
   "HEALPix: A Framework for High-Resolution Discretization and Fast Analysis of Data Distributed on the Sphere",
   *The Astrophysical Journal*, **622**, 759
   
   DOI: `10.1086/427976 <https://doi.org/10.1086/427976>`_
   
   ArXiv: `astro-ph/0409513 <https://arxiv.org/abs/astro-ph/0409513>`_

Software and Tools
------------------

.. [HEALPix] HEALPix software package
   
   Website: `https://healpix.jpl.nasa.gov/ <https://healpix.jpl.nasa.gov/>`_
   
   GitHub: `https://github.com/healpy/healpy <https://github.com/healpy/healpy>`_

.. [NPIPE] Planck NPIPE pipeline
   
   GitHub: `https://github.com/planck-npipe/toast-npipe <https://github.com/planck-npipe/toast-npipe>`_

.. [TOAST] Time Ordered Astrophysics Scalable Tools (TOAST)
   
   GitHub: `https://github.com/hpc4cmb/toast <https://github.com/hpc4cmb/toast>`_
   
   Documentation: `https://toast-cmb.readthedocs.io/ <https://toast-cmb.readthedocs.io/>`_

Data Archives
-------------

Planck Legacy Archive (PLA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Official repository for Planck mission data products:

* Website: `https://pla.esac.esa.int/ <https://pla.esac.esa.int/>`_
* Contains: Maps, power spectra, RIMO files, beam models, simulations
* Data releases: PR1 (2013), PR2 (2015), PR3 (2018), PR4/NPIPE (2020)

NASA Legacy Archive for Microwave Background Data Analysis (LAMBDA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternative archive with mirror of Planck data:

* Website: `https://lambda.gsfc.nasa.gov/ <https://lambda.gsfc.nasa.gov/>`_
* Includes: Planck, WMAP, COBE, and other CMB experiments

European Space Agency (ESA) Science Archives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Planck mission page: `https://www.cosmos.esa.int/web/planck <https://www.cosmos.esa.int/web/planck>`_
* Technical documentation, mission overview, publications

Related Software Packages
--------------------------

CMB Power Spectrum Estimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **PolSpice**: Pseudo-Cℓ estimator for polarized maps
  
  Website: `http://www2.iap.fr/users/hivon/software/PolSpice/ <http://www2.iap.fr/users/hivon/software/PolSpice/>`_

* **NaMaster**: Pseudo-Cℓ estimator with mode coupling matrices
  
  GitHub: `https://github.com/LSSTDESC/NaMaster <https://github.com/LSSTDESC/NaMaster>`_

* **MASTER**: Original pseudo-Cℓ method implementation
  
  Reference: Hivon et al. 2002, ApJ, 567, 2

CMB Map Making and Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Planck Level S**: Simulation and analysis pipeline
  
  Code available through PLA

* **Commander**: Global CMB analysis framework
  
  Used for Planck component separation

* **SMICA**: Spectral Matching Independent Component Analysis
  
  Planck CMB map reconstruction method

Python Scientific Stack
^^^^^^^^^^^^^^^^^^^^^^^^

Core dependencies used by qp-planck:

* **NumPy**: `https://numpy.org/ <https://numpy.org/>`_
* **SciPy**: `https://scipy.org/ <https://scipy.org/>`_
* **Astropy**: `https://www.astropy.org/ <https://www.astropy.org/>`_
* **matplotlib**: `https://matplotlib.org/ <https://matplotlib.org/>`_
* **mpi4py**: `https://mpi4py.readthedocs.io/ <https://mpi4py.readthedocs.io/>`_

Documentation and Tutorials
----------------------------

* **Planck Explanatory Supplement**: Comprehensive technical documentation
  
  `https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Explanatory_Supplement <https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Explanatory_Supplement>`_

* **HEALPix Primer**: Introduction to HEALPix pixelization
  
  Included with HEALPix software package

* **CMB-S4 Science Book**: Future CMB experiments and analysis methods
  
  ArXiv: `1907.04473 <https://arxiv.org/abs/1907.04473>`_

Citing qp-planck
----------------

If you use **qp-planck** in your research, please cite:

1. The QuickPol paper [Hivon2017]_
2. The appropriate Planck beam paper(s) [Planck2016_IV]_ and/or [Planck2016_VII]_
3. The NPIPE data release paper [Planck2020_LVII]_ if using NPIPE data
4. This software (include GitHub URL and version number)

Example BibTeX entries:

.. code-block:: bibtex

   @ARTICLE{Hivon2017,
       author = {{Hivon}, E. and {Mottet}, S. and {Ponthieu}, N.},
       title = "{QuickPol: Fast calculation of effective beam matrices for CMB polarization}",
       journal = {A\&A},
       year = 2017,
       volume = 598,
       eid = {A25},
       pages = {A25},
       doi = {10.1051/0004-6361/201629626},
       archivePrefix = {arXiv},
       eprint = {1606.09313},
   }
   
   @ARTICLE{Planck2020_LVII,
       author = {{Planck Collaboration}},
       title = "{Planck intermediate results. LVII. Joint Planck LFI and HFI data processing}",
       journal = {A\&A},
       year = 2020,
       volume = 643,
       eid = {A42},
       pages = {A42},
       doi = {10.1051/0004-6361/202038073},
   }
   
   @SOFTWARE{qp_planck,
       author = {{qp-planck Contributors}},
       title = "{qp-planck: QuickPol beam windows for Planck NPIPE}",
       howpublished = {\\url{https://github.com/paganol/qp_planck}},
       year = 2025,
   }

Contact and Support
-------------------

* **GitHub Issues**: `https://github.com/paganol/qp_planck/issues <https://github.com/paganol/qp_planck/issues>`_
* **Planck Help Desk**: For questions about Planck data products
* **HEALPix Forum**: For HEALPix-related questions

Acknowledgments
---------------

This software is adapted from the Planck NPIPE pipeline. We acknowledge:

* The Planck Collaboration for mission data and original analysis tools
* ESA for Planck mission support
* The HEALPix team for spherical harmonic software
* Contributors to NumPy, SciPy, Astropy, and other open-source tools

The Planck mission was supported by ESA, NASA, and numerous national space agencies and research institutions.
