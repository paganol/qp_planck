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
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

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
    blmfile: str = "blm_{}.fits",
    outdir: str = "../quickpol_output",
    indir: Optional[str] = None,
    smax: int = 6,
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
    do_plot: bool = False,
    # Optional YAML configuration
    config: Optional[Union[str, Path, Mapping]] = None,
) -> None:
    """
    Run the full QuickPol-style pipeline:
      hitmaps + blm  →  beam matrices (NPZ)  →  FITS B_ell / W_ell.

    You can configure it either:
      - by passing keyword arguments directly, or
      - by providing a YAML file via `config`, whose values will be used
        as defaults and overridden by explicit keyword arguments.

    Parameters
    ----------
    detpairs : sequence of (str, str)
        List of detector or detset pairs, e.g. [('143A', '143A'), ('143A', '143B')].
    rimo, rimo_lfi, rimo_hfi : str, optional
        Paths to RIMO FITS files. You can pass a single merged RIMO via `rimo`,
        or separate ones via `rimo_lfi` and/or `rimo_hfi`. These may also be
        supplied in the YAML config:
          rimo:     /path/to/RIMO.fits
          rimo_lfi: /path/to/RIMO_LFI_...
          rimo_hfi: /path/to/RIMO_HFI_...
    blmfile : str, optional
        Template for beam multipole files, e.g. 'blm_{}.fits'.
    outdir : str, optional
        Output directory for NPZ and FITS files.
    indir : str, optional
        Input directory for NPZ files when calling mat2fits.
        If None, defaults to `outdir`.
    smax : int, optional
        Maximum spin used in the QuickPol computation.
    spin_ref : {'Dxx', 'Pxx'}, optional
        Reference frame for the spin moments (as in original code).
    blm_ref : {'Dxx', 'Pxx'}, optional
        Reference frame for the input blm.
    angle_shift : float, optional
        Extra angle shift (degrees) applied in the beam matrix step.
    force_det : str, optional
        Force a particular detector label in filenames (QuickPol convention).
    release : str, optional
        Release tag used in output filenames (e.g. 'npipe6v20').
    rhobeam : {'IMO', 'Ideal'} or float, optional
        Beam cross-polarization model flag (passed to hmap2mat).
    rhohit : {'IMO', 'Ideal'} or float, optional
        Hit cross-polarization model flag (passed to hmap2mat).
    test : bool, optional
        If True, only sample a small fraction of pixels (debug / test mode).
    planet : str, optional
        Planet name used in deconvolution (if any, e.g. 'Saturn').
    conserve_memory : bool, optional
        If True, process maps in smaller chunks.
    overwrite : bool, optional
        Overwrite existing NPZ / FITS files.
    blfile, blTEBfile, wlfile : bool, optional
        Control which FITS outputs are produced (B_ell(T), B_ell(TEB), W_ell).
    do_plot : bool, optional
        If True, also produce PNG plots of the window functions.
    config : str or Path or Mapping, optional
        YAML file path or dict providing defaults for all the above.
        Explicit keyword arguments always take precedence.

    Notes
    -----
    YAML example:

        rimo_lfi: /path/to/RIMO_LFI_npipe5_symmetrized.fits
        rimo_hfi: /path/to/RIMO_HFI_npipe5v16_symmetrized.fits
        outdir: ./quickpol_output
        smax: 6
        release: npipe6v20
        rhobeam: IMO
        rhohit: IMO
        detpairs:
          - [143GHz, 143GHz]
          - [143A, 143A]
          - [143A, 143B]
    """
    # ------------------------------------------------------------------
    # Merge YAML config (if any) with explicit kwargs
    # ------------------------------------------------------------------
    cfg: Dict = {}
    if config is not None:
        cfg = _load_yaml_config(config)

    # Explicit kwargs override config values
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
    outdir = cfg.get("outdir", outdir)
    indir = cfg.get("indir", indir if indir is not None else outdir)
    smax = int(cfg.get("smax", smax))
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
    do_plot = bool(cfg.get("do_plot", do_plot))

    # detpairs can also come from config if not passed
    if not detpairs and "detpairs" in cfg:
        # YAML will give list of lists or tuples
        detpairs = [tuple(p) for p in cfg["detpairs"]]

    if not detpairs:
        raise ValueError("No detector pairs specified (detpairs is empty).")

    # ------------------------------------------------------------------
    # Main loop over detector pairs
    # ------------------------------------------------------------------
    for ipair, detpair in enumerate(detpairs):
        if ipair % ntask != rank:
            continue

        detset1, detset2 = detpair
        print(prefix, "Processing pair:", detset1, "x", detset2, flush=True)

        # 1) Build beam matrix NPZ via hmap2mat (former program1)
        hmap2mat(
            rimo_dict,
            detpair,
            blmfile,
            outdir,
            smax,
            spin_ref,
            blm_ref,
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

        # 2) Convert NPZ → FITS window functions via mat2fits (former q2f)
        mat2fits(
            indir,
            outdir,
            detpair,
            smax,
            release=release,
            full=not test,
            blfile=blfile,
            blTEBfile=blTEBfile,
            wlfile=wlfile,
            overwrite=overwrite,
            do_plot=do_plot,
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
    # detpairs will be taken from config inside run_qp_pipeline
    run_qp_pipeline(detpairs=(), config=config)


__all__ = ["run_qp_pipeline", "run_qp_from_yaml"]