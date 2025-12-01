"""
qp_planck
=========

DISCLAIMER
----------
This package contains code adapted from the Planck NPIPE pipeline:
    https://github.com/planck-npipe/toast-npipe/tree/master

The original implementation relied on TOAST and internal Planck tooling.
This version removes TOAST dependencies and restructures the code for
lightweight, standalone use in beam window-function calculations.

Description
-----------
Tools for computing beam window functions for Planck NPIPE using
QuickPol-inspired logic, and for exporting those results to FITS
B_ell / W_ell window function files.

Public API:
- load_RIMO()             → load instrument model (RIMO) FITS files
- list_planck()           → list detectors and detector sets
- detector_weights        → dictionary of detector white-noise weights
- qp_file()               → helper for producing QuickPol-style filenames
- qp                      → QuickPol beam matrix computation driver
- qp2fits                 → conversion of QuickPol NPZ outputs to FITS windows

Modules:
- qp.py        : main QuickPol driver (beam matrices, beam_mat NPZ)
- qp2fits.py   : NPZ → FITS B_ell / W_ell converter
- utilities.py : RIMO loading, detector lists, filename helpers
"""

from .utilities import (
    load_RIMO,
    list_planck,
    detector_weights,
    qp_file,
)

# Expose the qp and qp2fits modules so that users can do:
#   from qp_planck import qp, qp2fits
from . import qp
from . import qp2fits

__all__ = [
    "load_RIMO",
    "list_planck",
    "detector_weights",
    "qp_file",
    "qp",
    "qp2fits",
]
