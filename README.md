# qp_planck

QuickPol-based beam window function computation for Planck NPIPE.

This package wraps a lightly adapted version of the Planck NPIPE QuickPol
pipeline to:

- build **effective beam matrices** for Planck detector / detset combinations;
- read **RIMO** (reduced instrument model) files for LFI and HFI;
- read or generate **beam multipoles** (Bₗₘ);
- compress the result into **NPZ beam matrix files**;
- convert those NPZ files into **FITS Bℓ / Wℓ window functions**.

---

## Disclaimer

This code is **adapted from the Planck NPIPE pipeline**:

https://github.com/planck-npipe/toast-npipe/tree/master

The original implementation relied on TOAST and internal Planck tooling.
Here, the code has been minimally modified to remove TOAST dependencies and to
run as a standalone Python package.

---

## Installation

You can install `qp_planck` locally for development:

```bash
git clone https://github.com/paganol/qp_planck.git
cd qp_planck
pip install -e .
```

This will install all Python modules in editable mode.  
Make sure you have the following dependencies available:

- numpy  
- scipy  
- healpy  
- astropy  
- pyyaml  
- mpi4py (optional, only if using MPI execution)

---

## Package layout

After installation, the main package is:

```
qp_planck/
    qp_hmap2mat.py    # QuickPol driver: builds beam matrices and NPZ files (map → matrix)
    qp_mat2fits.py    # NPZ → FITS (Bℓ / Wℓ) window-function converter
    qp_pipeline.py    # High-level wrapper: runs hmap2mat + mat2fits; supports YAML configs
    utilities.py      # RIMO loading, detector lists, filename helpers, weights
    __init__.py       # Public API: load_RIMO, list_planck, qp_file,
                      #             hmap2mat, mat2fits, run_qp_pipeline
scripts/
    example_planck.yaml   # Example YAML configuration for running the full pipeline
    make_detpairs.py      # Helper for generating Planck detset / detector pairs
```

---

## Running modes

You can run the pipeline in **three** distinct ways.

### 1. Run using a YAML configuration (recommended)

You can drive the full pipeline (maps → matrices → FITS windows) with one command:

```bash
python -m qp_planck.qp_pipeline scripts/example_planck.yaml
```

The YAML file controls:

- list of detector pairs  
- path to RIMO files  
- blm file pattern  
- output directory  
- smax, lmax, release tag  
- whether to run hmap2mat, mat2fits, or both  

This is the cleanest way to reproduce full Planck-level outputs.

---

### 2. Run directly from Python

```python
from qp_planck import run_qp_pipeline

run_qp_pipeline(
    detpairs=[("100GHz", "100GHz")],
    rimo_lfi="path/to/RIMO_LFI.fits",
    rimo_hfi="path/to/RIMO_HFI.fits",
    outdir="output/",
    run_hmap2mat=True,
    run_mat2fits=True
)
```

---

### 3. CLI: run low-level steps manually

#### Generate NPZ beam matrices:

```bash
python -m qp_planck.qp_hmap2mat
```

#### Convert NPZ → FITS:

```bash
python -m qp_planck.qp_mat2fits
```

These follow the original Planck QuickPol tools, with file naming via `qp_file`.

---

## References

If you use this package, please cite:

**QuickPol**  
E. Hivon, S. Mottet, N. Ponthieu,  
*QuickPol: Fast calculation of effective beam matrices for CMB polarization*,  
A&A 598, A25 (2017).

**Planck NPIPE**  
Planck Collaboration Int. LVII, *NPIPE data processing pipeline*,
A&A 643, A42 (2020).
