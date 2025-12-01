# qp_planck

QuickPol-based beam window function computation for Planck NPIPE.

This package wraps a lightly adapted version of the Planck NPIPE QuickPol
pipeline to:

- build **effective beam matrices** for Planck detector / detset combinations;
- read **RIMO** (reduced instrument model) files for LFI and HFI;
- read or generate **beam multipoles** \(B\_{ℓm}\);
- compress the result into **NPZ beam matrix files**;
- convert those NPZ files into **FITS B\_ℓ / W\_ℓ window functions**.

---

## Disclaimer

This code is **adapted from the Planck NPIPE pipeline**:

> https://github.com/planck-npipe/toast-npipe/tree/master

The original implementation relied on TOAST and internal Planck tooling.
Here, the code has been minimally modified to remove TOAST dependencies and to
run as a standalone Python package.

If you use this in a scientific context, please also cite the QuickPol paper:

> E. Hivon, S. Mottet, N. Ponthieu,  
> *QuickPol: Fast calculation of effective beam matrices for CMB polarization*,  
> A&A 598, A25 (2017)

---

## Package layout

After installation, the main package is:

```text
qp_planck/
    qp.py        # QuickPol driver: builds beam matrices and NPZ files
    qp2fits.py   # NPZ → FITS Bℓ / Wℓ converter
    utilities.py # RIMO loading, detector lists, filenames, weights
    __init__.py  # public API: load_RIMO, list_planck, qp_file, qp, qp2fits, ...