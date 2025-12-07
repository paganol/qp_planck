"""
qp_pipeline.py
==============

High-level QuickPol pipeline wrapper for Planck NPIPE.

This module provides a convenience wrapper around:

  - qp_hmap2mat.hmap2mat : build effective beam matrices (NPZ)
  - qp_mat2fits.mat2fits : convert those matrices to FITS B_ell / W_ell

It also supports loading all options from a YAML configuration file.

DISCLAIMER
----------
This file contains code adapted from the Planck NPIPE pipeline:
    https://github.com/planck-npipe/toast-npipe/tree/master

The original implementation relied on TOAST and internal Planck tooling.
This version removes TOAST dependencies and is intended to be used as
a standalone step for beam window-function calculations.
"""

# Optional MPI via mpi4py; fall back to serial execution if not available.
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    ntask = comm.size
    rank = comm.rank
except Exception:  # pragma: no cover - serial fallback
    comm = None
    ntask = 1
    rank = 0

prefix = f"{rank:04d} :"

import os
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Tuple, Union

try:
    import yaml
except Exception:  # pragma: no cover - optional dependency
    yaml = None  # type: ignore

from .utilities import load_RIMO
from .qp_hmap2mat import hmap2mat
from .qp_mat2fits import mat2fits


# -----------------------------------------------------------------------------#
# YAML helpers                                                                 #
# -----------------------------------------------------------------------------#


def _load_yaml_config(config: Union[str, Path, Mapping]) -> Dict:
    """
    Load a configuration dictionary from a YAML file or mapping.

    Parameters
    ----------
    config : str or Path or Mapping
        If str / Path, path to a YAML file.
        If a mapping (dict-like), it is returned as a plain dict.

    Returns
    -------
    dict
        Configuration dictionary.
    """
    if isinstance(config, Mapping):
        return dict(config)

    path = Path(config)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    if yaml is None:
        raise ImportError(
            "PyYAML is required to read a YAML config file, "
            "but could not be imported."
        )

    with path.open("r") as f:
        cfg = yaml.safe_load(f) or {}
    if not isinstance(cfg, dict):
        raise ValueError(f"YAML config must produce a mapping, got {type(cfg)}")
    return cfg


# -----------------------------------------------------------------------------#
# RIMO helper                                                                  #
# -----------------------------------------------------------------------------#


def _load_rimo_from_config(cfg: Mapping) -> Dict:
    """
    Load and merge LFI/HFI RIMO files based on config entries.

    Expected keys (any combination):
      - rimo_lfi : path to LFI RIMO FITS
      - rimo_hfi : path to HFI RIMO FITS
      - rimo     : single RIMO FITS (if you have a merged file already)

    Returns
    -------
    dict
        Merged RIMO dictionary (keys are detector names).
    """
    rimo = {}

    rimo_single = cfg.get("rimo")
    rimo_lfi = cfg.get("rimo_lfi")
    rimo_hfi = cfg.get("rimo_hfi")

    if rimo_single is not None:
        print(prefix, "Loading RIMO from", rimo_single, flush=True)
        rimo.update(load_RIMO(rimo_single))

    if rimo_lfi is not None:
        print(prefix, "Loading LFI RIMO from", rimo_lfi, flush=True)
        rimo.update(load_RIMO(rimo_lfi))

    if rimo_hfi is not None:
        print(prefix, "Loading HFI RIMO from", rimo_hfi, flush=True)
        rimo.update(load_RIMO(rimo_hfi))

    if not rimo:
        raise ValueError(
            "No RIMO file specified. Provide 'rimo', or 'rimo_lfi'/'rimo_hfi' "
            "in the config or function arguments."
        )

    return rimo


# -----------------------------------------------------------------------------#
# Main high-level pipeline                                                     #
# -----------------------------------------------------------------------------#


def run_qp_pipeline(
    detpairs: Sequence[Tuple[str, str]],
    *,
    # RIMO files (can also be provided via YAML)
    rimo: Optional[str] = None,
    rimo_lfi: Optional[str] = None,
    rimo_hfi: Optional[str] = None,
    # Beam & hitmap configuration
    blmfile: str = "inputs/beams/blm_{}.fits",
    momentsdir: str = "inputs/polmoments",
    outdir: str = "../quickpol_output",
    indir: Optional[str] = None,
    smax: int = 6,
    nside: Optional[str] = None,
    lmax: Optional[str] = None,
    mmax: int = 10,
    spin_ref: str = "Pxx",
    blm_ref: str = "Dxx",
    angle_shift: float = 0.0,
    force_det: Optional[str] = None,
    release: str = "npipe6v20",
    rhobeam: str = "IMO",
    rhohit: str = "IMO",
    test: bool = False,
    planet: str = "",
    conserve_memory: bool = True,
    overwrite: bool = False,
    # FITS / window-function options
    blfile: bool = True,
    blTEBfile: bool = True,
    wlfile: bool = True,
    # Optional YAML configuration
    config: Optional[Union[str, Path, Mapping]] = None,
) -> None:
    """
    Run the full QuickPol pipeline for a list of detector pairs.

    This function coordinates all major QuickPol steps:

        hitmaps + b_lm  →  beam matrices (NPZ)  →  B_ell / W_ell FITS products

    It first loads configuration parameters (from explicit keyword arguments
    or an optional YAML file), then loops over all detector–set pairs and:

      1. Calls ``hmap2mat`` to compute the beam matrix NPZ file.
      2. Calls ``mat2fits`` to generate FITS window-function products
         (B_ell, B_ell^TEB, W_ell) and optional diagnostic plots.

    Keyword arguments always override the YAML configuration.

    Parameters
    ----------
    detpairs : sequence of (str, str)
        List of detector or detector-set pairs, e.g.
        ``[("143-1a", "143-1b"), ("143-2a", "143-2a")]``.
        If empty, the function attempts to load ``detpairs`` from the YAML config.

    config : str or Path or mapping, optional
        YAML configuration file or a Python dictionary.  Values in the YAML file
        serve as defaults and are overridden by explicit keyword arguments.

    rimo, rimo_lfi, rimo_hfi : str, optional
        Paths to RIMO files.  You may supply one merged RIMO via ``rimo`` or
        separate LFI/HFI RIMOs.  These may also be specified in the YAML file.

    blmfile : str, optional
        Template for beam multipole filenames (e.g. ``"blm_{}.fits"``).

    momentsdir : str, optional
        Directory containing detector hit-moment maps (polmoments files).

    outdir : str, optional
        Directory where NPZ and FITS outputs are written.

    indir : str, optional
        Directory where NPZ files are read by ``mat2fits``.
        If omitted, defaults to ``outdir``.

    smax : int, optional
        Maximum spin in the hit matrix and beam computation.

    nside : int or None, optional
        HEALPix resolution used when building hit matrices.  If None, each
        detector pair determines its own ``nside`` through QuickPol conventions.

    lmax : int or None, optional
        Maximum multipole.  If None, defaults to ``4 * nside`` per detector pair.

    mmax : int, optional
        Maximum |m| index to load from the beam b_lm files.

    spin_ref : str, optional
        QuickPol reference spin convention (e.g. ``"Pxx"`` or ``"Dxx"``).

    blm_ref : str, optional
        Reference b_lm set inside the beam file.

    angle_shift : float, optional
        Additional polarization angle shift (degrees).

    force_det : str, optional
        Override detector name when creating output filenames.

    release : str, optional
        Pipeline version / release tag encoded in filenames.

    rhobeam, rhohit : {"IMO", "Ideal"} or float, optional
        Cross-polarization model parameters.

    test : bool, optional
        If True, run in reduced-pixel test mode.

    planet : str, optional
        Planet label used in the beam matrix computation.

    conserve_memory : bool, optional
        Reduce memory usage at the cost of extra CPU.

    overwrite : bool, optional
        If True, overwrite existing NPZ and FITS files.

    blfile, blTEBfile, wlfile : bool, optional
        Flags controlling which FITS products are written.

    Notes
    -----
    • Argument resolution order is:

        explicit keyword arguments  >  YAML config values  >  built-in defaults

    • Parallel execution: each MPI rank processes a subset of detector pairs.

    • This function does not return anything: all results are written to disk.
    """

    # ------------------------------------------------------------------
    # Merge YAML config (if any) with explicit kwargs
    # ------------------------------------------------------------------
    cfg: Dict = {}
    if config is not None:
        cfg = _load_yaml_config(config)

    # Helper: explicit kwargs override config values
    def _get(name, current):
        return current if current is not None else cfg.get(name)

    # Resolve RIMO paths from config + args
    rimo_paths_cfg: Dict[str, Optional[str]] = {
        "rimo": _get("rimo", rimo),
        "rimo_lfi": _get("rimo_lfi", rimo_lfi),
        "rimo_hfi": _get("rimo_hfi", rimo_hfi),
    }
    # drop None values
    rimo_paths_cfg = {k: v for k, v in rimo_paths_cfg.items() if v is not None}

    if not rimo_paths_cfg:
        # maybe config only
        rimo_dict = _load_rimo_from_config(cfg)
    else:
        # build a small cfg-like dict and reuse helper
        rimo_dict = _load_rimo_from_config(rimo_paths_cfg)

    # other options possibly overridden by config
    blmfile = cfg.get("blmfile", blmfile)
    momentsdir = cfg.get("momentsdir", momentsdir)
    outdir = cfg.get("outdir", outdir)
    indir = cfg.get("indir", indir if indir is not None else outdir)
    smax = int(cfg.get("smax", smax))
    nside = cfg.get("nside", nside)
    nside = int(nside) if nside is not None else None
    lmax = cfg.get("lmax", lmax)
    lmax = int(lmax) if lmax is not None else None
    spin_ref = cfg.get("spin_ref", spin_ref)
    blm_ref = cfg.get("blm_ref", blm_ref)
    angle_shift = float(cfg.get("angle_shift", angle_shift))
    force_det = cfg.get("force_det", force_det)
    release = cfg.get("release", release)
    rhobeam = cfg.get("rhobeam", rhobeam)
    rhohit = cfg.get("rhohit", rhohit)
    test = bool(cfg.get("test", test))
    planet = cfg.get("planet", planet)
    conserve_memory = bool(cfg.get("conserve_memory", conserve_memory))
    overwrite = bool(cfg.get("overwrite", overwrite))
    blfile = bool(cfg.get("blfile", blfile))
    blTEBfile = bool(cfg.get("blTEBfile", blTEBfile))
    wlfile = bool(cfg.get("wlfile", wlfile))

    # detpairs can also come from config if not passed
    if (not detpairs) and ("detpairs" in cfg):
        detpairs = [tuple(p) for p in cfg["detpairs"]]

    if not detpairs:
        raise ValueError("No detector pairs specified (detpairs is empty).")

    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)

    # ------------------------------------------------------------------
    # Main loop over detector pairs
    # ------------------------------------------------------------------
    for ipair, detpair in enumerate(detpairs):
        if ipair % ntask != rank:
            continue

        detset1, detset2 = detpair
        print(prefix, "Processing pair:", detset1, "x", detset2, flush=True)

        # 1) Build beam matrix NPZ via hmap2mat
        hmap2mat(
            rimo_dict,
            detpair,
            blmfile,
            momentsdir,
            outdir,
            smax,
            spin_ref,
            blm_ref,
            nside=nside,
            lmax=lmax,
            mmax=mmax,
            angle_shift=angle_shift,
            force_det=force_det,
            release=release,
            rhobeam=rhobeam,
            rhohit=rhohit,
            test=test,
            planet=planet,
            conserve_memory=conserve_memory,
            overwrite=overwrite,
        )

        # 2) Convert NPZ → FITS window functions via mat2fits
        mat2fits(
            indir,
            outdir,
            detpair,
            smax,
            lmax=lmax,
            release=release,
            full=not test,
            blfile=blfile,
            blTEBfile=blTEBfile,
            wlfile=wlfile,
            overwrite=overwrite,
        )


# -----------------------------------------------------------------------------#
# Convenience function: run from YAML only                                     #
# -----------------------------------------------------------------------------#


def run_qp_from_yaml(config: Union[str, Path, Mapping]) -> None:
    """
    Run the QuickPol pipeline using only a YAML (or dict) configuration.

    This is just a thin wrapper around `run_qp_pipeline(config=...)`.

    Parameters
    ----------
    config : str or Path or Mapping
        YAML file path or dict containing all options, including
        a 'detpairs' entry.
    """
    run_qp_pipeline(detpairs=(), config=config)


__all__ = ["run_qp_pipeline", "run_qp_from_yaml"]


# -----------------------------------------------------------------------------#
# CLI entry point: `python -m qp_planck.qp_pipeline <config.yaml>`             #
# -----------------------------------------------------------------------------#

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python -m qp_planck.qp_pipeline <config.yaml>")
        raise SystemExit(1)

    config_path = sys.argv[1]
    run_qp_from_yaml(config_path)