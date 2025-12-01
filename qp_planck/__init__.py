"""
qp_planck
=========

DISCLAIMER
----------
This package contains code **adapted from the Planck NPIPE pipeline**:
    https://github.com/planck-npipe/toast-npipe/tree/master

The original implementation relied on TOAST and internal Planck tooling.
This version removes TOAST dependencies and restructures the code for
lightweight, standalone use in beam window-function calculations.

Description
-----------
Tools for computing beam window functions for Planck NPIPE using
QuickPol-inspired logic.

Public API:
- load_RIMO()    → load instrument model (RIMO) FITS files
- list_planck()  → list detectors and detector sets
- detector_weights
                 → dictionary of detector white-noise weights
- qp_file()      → helper for producing QuickPol-style filenames
- hmap2mat()     → build beam matrices from hitmaps / blm inputs
- mat2fits()     → convert beam matrices (NPZ) to B_ell / W_ell FITS files

Modules:
- utilities.py     : RIMO loading, detector lists, filename helpers
- qp_hmap2mat.py   : main QuickPol-style driver (hitmap → matrix)
- qp_mat2fits.py   : matrix → FITS window-function converter
"""

from .utilities import (
    load_RIMO,
    list_planck,
    detector_weights,
    qp_file,
)

from .qp_hmap2mat import hmap2mat
from .qp_mat2fits import mat2fits

__all__ = [
    "load_RIMO",
    "list_planck",
    "detector_weights",
    "qp_file",
    "hmap2mat",
    "mat2fits",
]