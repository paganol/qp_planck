Physics Background
==================

This page provides the theoretical foundation for the QuickPol beam window
function formalism implemented in **qp-planck**.

Beam Window Functions
----------------------

In CMB observations, the measured angular power spectrum :math:`C_\ell^{\text{obs}}`
is related to the true sky power spectrum :math:`C_\ell^{\text{sky}}` through
the beam window function :math:`W_\ell`:

.. math::

   C_\ell^{\text{obs}} = C_\ell^{\text{sky}} \cdot W_\ell

For polarization-sensitive observations, this generalizes to a **beam transfer matrix**
that couples temperature and polarization modes:

.. math::

   \begin{pmatrix}
   C_\ell^{TT,\text{obs}} \\
   C_\ell^{EE,\text{obs}} \\
   C_\ell^{BB,\text{obs}} \\
   C_\ell^{TE,\text{obs}}
   \end{pmatrix}
   =
   \mathbf{W}_\ell
   \begin{pmatrix}
   C_\ell^{TT,\text{sky}} \\
   C_\ell^{EE,\text{sky}} \\
   C_\ell^{BB,\text{sky}} \\
   C_\ell^{TE,\text{sky}}
   \end{pmatrix}

The QuickPol Method
-------------------

QuickPol (Hivon, Mottet & Ponthieu 2017) is a fast algorithm for computing
these beam matrices accounting for:

* **Asymmetric beams** — Non-circular beam shapes
* **Detector mismatch** — Different beams for different detectors
* **Polarization cross-talk** — Leakage between T, E, B modes
* **Scanning strategy** — Full integration over the mission timeline

The key innovation is representing the beam convolution in harmonic space
using **spin-weighted spherical harmonics**.

Spin-Weighted Harmonics
------------------------

A polarized beam is described by beam multipoles :math:`b_{\ell m}^{(s)}` for
each spin component :math:`s`:

* :math:`s = 0`: Temperature beam
* :math:`s = \pm 2`: Polarization beams (E/B)

The effective beam matrix element for multipole :math:`\ell` and input power
spectrum type :math:`XY` is:

.. math::

   W_\ell^{XY} = \sum_{d_1, d_2} \sum_{s=-s_{\max}}^{s_{\max}} 
   w_{d_1} w_{d_2} \, h_{d_1,d_2}^{(s)} \,
   B_{\ell,s}^{d_1} \, (B_{\ell,s}^{d_2})^*

where:

* :math:`d_1, d_2` are detectors in the map
* :math:`w_d` are detector weights (e.g., inverse noise variance)
* :math:`h_{d_1,d_2}^{(s)}` are **hit matrix** elements from the scanning pattern
* :math:`B_{\ell,s}^d` are beam matrix elements constructed from :math:`b_{\ell m}^d`

Hit Matrix Computation
-----------------------

The hit matrix :math:`h_{d_1,d_2}^{(s)}` encodes how the scanning strategy
couples different spin components:

.. math::

   h_{d_1,d_2}^{(s)} = \sum_{t \in \text{observations}} 
   e^{i s (\psi_{t,d_1} - \psi_{t,d_2})}

where :math:`\psi_{t,d}` is the polarization angle of detector :math:`d` at
time :math:`t`.

In practice, this is computed from **spin moment maps**:

.. math::

   \text{spinmap}_s(\mathbf{n}) = \sum_{t: \hat{\mathbf{n}}_t = \mathbf{n}} 
   e^{i s \psi_t}

stored in HEALPix format for each detector.

Beam Matrix Elements
---------------------

From beam multipoles :math:`b_{\ell m}`, we construct the beam matrix using
the relation:

.. math::

   B_{\ell,s}^{(d)} = 
   \begin{pmatrix}
   b_{\ell,s}^{(0,d)} & b_{\ell,s-2}^{(0,d)} & b_{\ell,s+2}^{(0,d)} \\
   b_{\ell,s}^{(+,d)} & b_{\ell,s-2}^{(+,d)} & b_{\ell,s+2}^{(+,d)} \\
   b_{\ell,s}^{(-,d)} & b_{\ell,s-2}^{(-,d)} & b_{\ell,s+2}^{(-,d)}
   \end{pmatrix}

where:

* :math:`b_{\ell,m}^{(0)}`: temperature beam component
* :math:`b_{\ell,m}^{(\pm)}`: polarization beam components (E ± iB)

The beam multipoles incorporate:

1. **Detector angle** :math:`\phi_d` → rotation in :math:`m` via :math:`e^{im\phi_d}`
2. **Polarization efficiency** :math:`\rho_d` → cross-polar response
3. **Beam FWHM** and ellipticity from optical modeling

Polarization Efficiency
------------------------

For Planck detectors, two models are available:

**Ideal Model** (``rhobeam='Ideal'``):

.. math::

   \rho = 
   \begin{cases}
   1 & \text{for PSBs (Polarization Sensitive Bolometers)} \\
   0 & \text{for SWBs (Spider-Web Bolometers)}
   \end{cases}

**IMO Model** (``rhobeam='IMO'``):

Uses measured polarization efficiency from the Instrument Model (RIMO):

.. math::

   \rho = \frac{1 - \epsilon}{1 + \epsilon}

where :math:`\epsilon` is the cross-polarization fraction.

Output Products
---------------

The pipeline produces two types of window functions:

1. **Scalar Beam Windows** (:math:`B_\ell`)

   Diagonal approximation for temperature or single-field analyses:
   
   .. math::
   
      C_\ell^{\text{obs}} = C_\ell^{\text{sky}} \cdot B_\ell^2

   Stored as 1D arrays in FITS files (``Bl_*.fits``).

2. **Matrix Beam Windows** (:math:`W_\ell`)

   Full :math:`3 \times 3` coupling matrix for T, E, B:
   
   .. math::
   
      W_\ell = 
      \begin{pmatrix}
      W_\ell^{TT \to TT} & W_\ell^{TE \to TT} & W_\ell^{TB \to TT} \\
      W_\ell^{TE \to EE} & W_\ell^{EE \to EE} & W_\ell^{EB \to EE} \\
      W_\ell^{TB \to BB} & W_\ell^{EB \to BB} & W_\ell^{BB \to BB}
      \end{pmatrix}

   Includes cross-talk terms (e.g., :math:`W_\ell^{TE \to TT}` is temperature
   contamination from TE input).
   
   Stored in binary FITS tables (``Wl_*.fits``) with 9 columns per :math:`\ell`.

Numerical Implementation
-------------------------

**Key Parameters:**

:math:`\ell_{\max}`
   Maximum multipole (typically 4096 for HFI, 2048 for LFI)

:math:`s_{\max}`
   Maximum spin moment (default: 6, sufficient for Planck)

:math:`m_{\max}`
   Maximum azimuthal mode in :math:`b_{\ell m}` (typically 10)

**Memory Optimization:**

The code uses block-wise inversion of the hit matrix and optional
``conserve_memory`` mode to reduce RAM usage for large detector sets.

**MPI Parallelization:**

Detector pairs are distributed across MPI ranks, allowing linear scaling
to 100+ cores on HPC systems.

Further Reading
---------------

* **QuickPol Paper**: Hivon, E., Mottet, S., & Ponthieu, N. 2017, A&A, 598, A25
  
  *"QuickPol: Fast calculation of effective beam matrices for CMB polarization"*
  
  https://doi.org/10.1051/0004-6361/201629626

* **Planck Beam Papers**:
  
  * Planck Collaboration IV. 2016, A&A, 594, A4 (LFI beams)
  * Planck Collaboration VII. 2016, A&A, 594, A7 (HFI beams)

* **NPIPE Data Release**: Planck Collaboration Int. LVII. 2020, A&A, 643, A42

* **HEALPix**: Górski, K. M., et al. 2005, ApJ, 622, 759

See Also
--------

* :doc:`api/hmap2mat` — Implementation details of the beam matrix calculation
* :doc:`examples` — Practical examples with real Planck data
* :doc:`references` — Complete bibliography
