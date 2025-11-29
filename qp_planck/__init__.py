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
- load_RIMO()            → load instrument model (RIMO) FITS files
- list_planck()          → list detectors and detector sets
- detector_weights        → dictionary of detector white-noise weights
- qp_file()               → helper for producing QuickPol-style filenames
- qp                      → full QuickPol wrapper and beam matrix generation

Modules:
- qp.py          : main QuickPol driver
- utilities.py   : RIMO loading, detector lists, filename helpers
"""
from .utilities import (
    load_RIMO,
    list_planck,
    detector_weights,
    qp_file,
)

# Expose the qp module itself so that users can do:
#    from qp_planck import qp
from . import qp

__all__ = [
    "load_RIMO",
    "list_planck",
    "detector_weights",
    "qp_file",
    "qp",
]
