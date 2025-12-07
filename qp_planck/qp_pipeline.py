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
    nside: Optional[int] = None,
    lmax: Optional[int] = None,
    mmax: int = 10,
    spin_ref: str = "Pxx",
    blm_ref: str = "Dxx",
    masks: Optional[str] = None,
    masks_names: Optional[str] = None,
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

        (hitmaps + b_lm + masks)  →  beam matrices (NPZ)  →  B_ell / W_ell FITS products

    The function first loads configuration parameters—either from explicit
    keyword arguments or from an optional YAML configuration file—then loops
    over all detector pairs and:

      1. Calls ``hmap2mat`` to compute the effective beam matrix and metadata.
      2. Calls ``mat2fits`` to generate FITS window-function products:
         - B_ell (temperature or temperature+polarization),
         - B_ell^{TEB},
         - W_ell (cut-sky window functions),
         as well as optional diagnostic information.

    Explicit keyword arguments always override values provided in the YAML file.

    Parameters
    ----------
    detpairs : sequence of (str, str)
        List of detector pairs, e.g.
        ``[("143-1a", "143-1b"), ("217-5a", "217-5b")]``.
        If empty, ``detpairs`` may be loaded from the YAML configuration.

    config : str, Path, or mapping, optional
        YAML configuration file or a Python dictionary.  Settings from the YAML
        act as defaults and are overridden by explicit keyword arguments.

    rimo, rimo_lfi, rimo_hfi : str, optional
        Paths to RIMO (Radiometer Instrument Model) files.  Either one merged
        RIMO may be given via ``rimo`` or separate LFI/HFI RIMOs via
        ``rimo_lfi`` and ``rimo_hfi``.

    blmfile : str
        Template for beam multipole filenames, e.g. ``"blm_{}.fits"``.

    momentsdir : str
        Directory containing polarization-moment maps (polmoments) and hits.

    outdir : str
        Directory where `.npz` and FITS window-function files are written.

    indir : str or None, optional
        Directory from which ``mat2fits`` reads NPZ files.  
        If None, defaults to ``outdir``.

    smax : int
        Maximum spin used in hit/moment expansion.

    nside : int or None
        HEALPix resolution. If None, ``hmap2mat`` determines it per detector pair.

    lmax : int or None
        Maximum multipole. If None, defaults to ``4 * nside`` inside ``hmap2mat``.

    mmax : int, optional
        Maximum |m| index to load from each detector b_lm.

    spin_ref : str, optional
        QuickPol reference spin convention (e.g. ``"Pxx"`` or ``"Dxx"``).

    blm_ref : str, optional
        Reference beam inside the beam multipole files.

    masks : None, str, or dict, optional
        Mask specification:
          • ``None`` — full sky (no mask),  
          • ``str`` — filename of a mask to apply to all detector pairs,  
        Passed directly to ``hmap2mat`` and interpreted by ``get_all_masks``.

    masks_names : str or sequence of str, optional
        Optional human-readable mask names for metadata.

    angle_shift : float
        Additional rotation (degrees) applied to detector polarization angles.

    force_det : str or None
        Override detector name when generating filenames.

    release : str
        Processing or calibration release tag used for filenames.

    rhobeam, rhohit : {"IMO", "Ideal"} or float
        Cross-polarization or hit-matrix regularization parameters.

    test : bool
        If True, enable reduced test mode and skip expensive operations.

    planet : str
        Planet name used for beam normalization.

    conserve_memory : bool
        Aggressively free intermediate arrays inside ``hmap2mat`` and ``mat2fits``.

    overwrite : bool
        If True, overwrite existing NPZ or FITS outputs.  
        If False, existing files cause the corresponding steps to be skipped.

    blfile, blTEBfile, wlfile : bool
        Flags controlling which FITS window-function products are produced.

    Notes
    -----
    • Precedence of configuration values is:

          explicit keyword arguments  >  YAML config values  >  built-in defaults

    • Under MPI, each rank processes a subset of ``detpairs`` independently.

    • This function returns nothing: all results are written to disk.
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
    masks = cfg.get("masks", masks)    
    masks_names = cfg.get("masks_names", masks_names)
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
            masks=masks,
            masks_names=masks_names,
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