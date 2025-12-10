"""
Utility functions and data structures for Planck / QuickPol tools.

This module provides:

- RIMO (Reduced Instrument Model) loading and broadcasting.
- Convenience functions to list Planck detector sets.
- Detector weights for map combinations.
- A helper to construct QuickPol beam matrix filenames.

The code is written to be self-contained for use in a small qp_planck package.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np
from astropy.io import fits as pf

# ----------------------------------------------------------------------
# Basic constants and helpers
# ----------------------------------------------------------------------

#: Conversion factor from degrees to radians.
degree: float = np.pi / 180.0


class Timer:
    """
    Simple wall-clock timer utility.

    Used to roughly time operations and print a message when done.
    Compatible with the original usage pattern:

        timer = Timer()
        timer.start()
        ...
        timer.report_clear("message")
    """

    def __init__(self) -> None:
        self._t0: Optional[float] = None

    def start(self) -> None:
        """Start (or restart) the timer."""
        self._t0 = time.time()

    def report_clear(self, label: str) -> None:
        """
        Print the elapsed time since `start()` with a given label and clear the timer.

        Parameters
        ----------
        label : str
            Human-readable description of the measured operation.
        """
        if self._t0 is None:
            print(f"{label}: timer was not started")
            return
        dt = time.time() - self._t0
        print(f"{label}: {dt:.3f} s")
        self._t0 = None


@dataclass
class DetectorData:
    """
    Container for a single Planck detector's RIMO / focal plane information.

    Parameters
    ----------
    name : str
        Detector identifier (e.g. 'LFI27M', '143-1a', etc.).
    phi_uv : float
        Azimuth angle of the detector direction in degrees (UV frame).
    theta_uv : float
        Polar angle of the detector direction in degrees (UV frame).
    psi_uv : float
        Additional rotation angle in degrees (UV frame).
    psi_pol : float
        Polarization angle in degrees.
    epsilon : float
        Polarization efficiency (0–1).
    fsample : float
        Sample frequency in Hz.
    fknee : float
        1/f noise knee frequency in Hz.
    alpha : float
        1/f noise spectral index.
    net : float
        Noise Equivalent Temperature (NET).
    fwhm : float
        Beam full width at half maximum, typically in arcminutes.
    quat : np.ndarray
        Detector pointing quaternion (x, y, z, w) with scalar part last.
    """

    name: str
    phi_uv: float
    theta_uv: float
    psi_uv: float
    psi_pol: float
    epsilon: float
    fsample: float
    fknee: float
    alpha: float
    net: float
    fwhm: float
    quat: np.ndarray


# ----------------------------------------------------------------------
# nside and lmax utilities
# ----------------------------------------------------------------------

def detset2nside(detset):
    """Map a detset string to an appropriate nside."""
    if detset.startswith("0") or detset.startswith("LFI"):
        nside = 1024
    else:
        nside = 2048
    return nside

def detset2lmax(detset):
    """
    Map a detector set string to a default l_max.

    Parameters
    ----------
    detset : str
        Detector or detector set identifier.

    Returns
    -------
    int
        Default maximum multipole.
    """
    if detset.startswith("0") or detset.startswith("LFI"):
        lmax = 4 * 1024
    else:
        lmax = 4 * 2048
    return lmax


def detset2pol(detset):
    """
    Decide whether a detector set has polarization information.

    Parameters
    ----------
    detset : str

    Returns
    -------
    bool
        True if polarized, False otherwise.
    """
    if "545" in detset or "857" in detset or "LFI" in detset or "-" in detset:
        pol = False
    else:
        pol = True
    return pol


# ----------------------------------------------------------------------
# Quaternion utilities
# ----------------------------------------------------------------------


def quat_mult(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    """
    Multiply two quaternions with scalar part in the last component.

    Quaternion convention:
        q = [x, y, z, w]  with w the scalar part.

    Parameters
    ----------
    q1, q2 : np.ndarray
        Input quaternions, shape (4,).

    Returns
    -------
    np.ndarray
        Product quaternion q = q1 * q2, shape (4,).
    """
    x1, y1, z1, w1 = q1
    x2, y2, z2, w2 = q2

    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
    z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2

    return np.array([x, y, z, w])


class _QuaternionArrayNamespace:
    """
    Minimal stand-in for a 'quaternionarray' module used in some Planck code.

    Only provides a staticmethod `mult(q1, q2)` to match the old API.
    """

    @staticmethod
    def mult(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
        return quat_mult(q1, q2)


#: Namespace object to mimic 'qa.mult' used in legacy code.
qa = _QuaternionArrayNamespace()

#: Bore-sight rotation quaternion. Set to identity by default (no extra rotation).
SPINROT: np.ndarray = np.array([0.0, 0.0, 0.0, 1.0])


# ----------------------------------------------------------------------
# RIMO loading
# ----------------------------------------------------------------------


def load_RIMO(path: str, comm: Optional[Any] = None) -> Dict[str, DetectorData]:
    """
    Load and (optionally) broadcast the reduced instrument model (RIMO).

    The RIMO is the "focal plane database" for Planck: it describes
    detector positions, orientations, beam properties and noise parameters.

    Parameters
    ----------
    path : str
        Path to the FITS RIMO file.
    comm : MPI communicator, optional
        MPI communicator with attributes `rank`, and methods `Barrier()` and
        `bcast(obj, root=0)`. If provided, the RIMO dictionary is only
        loaded on rank 0 and then broadcast to all ranks.

    Returns
    -------
    dict
        Mapping from detector name (str) to a :class:`DetectorData` instance.
    """
    if comm is not None:
        comm.Barrier()

    timer = Timer()
    timer.start()

    RIMO: Dict[str, DetectorData] = {}
    is_root = comm is None or getattr(comm, "rank", 0) == 0

    if is_root:
        print(f"Loading RIMO from {path}", flush=True)
        hdulist = pf.open(path, "readonly")

        data = hdulist[1].data
        detectors = data.field("detector").ravel()
        phi_uvs = data.field("phi_uv").ravel()
        theta_uvs = data.field("theta_uv").ravel()
        psi_uvs = data.field("psi_uv").ravel()
        psi_pols = data.field("psi_pol").ravel()
        epsilons = data.field("epsilon").ravel()
        fsamples = data.field("f_samp").ravel()
        fknees = data.field("f_knee").ravel()
        alphas = data.field("alpha").ravel()
        nets = data.field("net").ravel()
        fwhms = data.field("fwhm").ravel()

        for i in range(len(detectors)):
            phi = phi_uvs[i] * degree
            theta = theta_uvs[i] * degree
            # Make sure we don't double-count psi rotation already included in phi
            psi = (psi_uvs[i] + psi_pols[i]) * degree - phi

            quat = np.zeros(4)
            # ZYZ convention conversion, scalar part at the end:
            # See e.g. http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf
            quat[3] = np.cos(0.5 * theta) * np.cos(0.5 * (phi + psi))  # scalar part
            quat[0] = -np.sin(0.5 * theta) * np.sin(0.5 * (phi - psi))
            quat[1] = np.sin(0.5 * theta) * np.cos(0.5 * (phi - psi))
            quat[2] = np.cos(0.5 * theta) * np.sin(0.5 * (phi + psi))

            # Apply the bore-sight rotation to the detector quaternion
            quat = qa.mult(SPINROT, quat)

            RIMO[detectors[i]] = DetectorData(
                name=detectors[i],
                phi_uv=phi_uvs[i],
                theta_uv=theta_uvs[i],
                psi_uv=psi_uvs[i],
                psi_pol=psi_pols[i],
                epsilon=epsilons[i],
                fsample=fsamples[i],
                fknee=fknees[i],
                alpha=alphas[i],
                net=nets[i],
                fwhm=fwhms[i],
                quat=quat,
            )

        hdulist.close()

    if comm is not None:
        # Broadcast to all ranks
        RIMO = comm.bcast(RIMO, root=0)

    if is_root:
        timer.report_clear("Load and broadcast RIMO")

    return RIMO


# ----------------------------------------------------------------------
# Planck detector lists
# ----------------------------------------------------------------------


def list_planck(
    detset: Union[int, str],
    good: bool = True,
    subset: Optional[int] = None,
    extend_857: bool = True,
    extend_545: bool = False,
) -> Union[List[str], int]:
    """
    Return lists of Planck detectors / horns for common detector sets.

    The logic follows the usual Planck conventions for LFI and HFI:

    Examples
    --------
    >>> list_planck(100)              # all 100 GHz detectors
    >>> list_planck("143A")           # subset 'A' at 143 GHz
    >>> list_planck("LFI")            # all LFI detectors
    >>> list_planck("PLANCK")         # all Planck detectors

    Special values
    --------------
    detset = "ROWS"
        Returns a list of lists that group detectors into rows, useful for
        plots or tables.
    detset = "LFI", "HFI", "PLANCK"
        Returns lists of all detectors in the respective instrument(s).

    Parameters
    ----------
    detset : int or str
        Detector set identifier, e.g. 30, "100", "143A", "LFI", "HFI", "PLANCK".
    good : bool, optional
        If True, only include "good" detectors where that distinction exists.
    subset : int, optional
        For some frequencies, can select specific subsets (1, 2, or 3)
        following Planck HFI conventions.
    extend_857 : bool, optional
        If True, at 857 GHz include extended set when `good` is False.
    extend_545 : bool, optional
        If True, at 545 GHz include extended set even if `good` is True.

    Returns
    -------
    list[str] or int
        A list of detector names, or `-1` if the detector set is unknown.
    """
    detectors: List[str] = []
    if subset is None:
        subset = 0

    # --- LFI channel selections ------------------------------------------------
    if detset in (30, "30", "030", "30GHz", "030GHz", "30A", "030A", "30B", "030B"):
        horns = range(27, 29)
        instrument = "LFI"
    elif detset in (44, "44", "044", "44GHz", "044GHz", "44A", "044A", "44B", "044B"):
        horns = range(24, 27)
        instrument = "LFI"
    elif detset in (70, "70", "070", "70GHz", "070GHz"):
        horns = range(18, 24)
        if subset == 1:
            horns = [18, 23]
        elif subset == 2:
            horns = [19, 22]
        elif subset == 3:
            horns = [20, 21]
        instrument = "LFI"
    elif detset in ["70A", "070A"]:
        horns = [18, 20, 23]
        instrument = "LFI"
    elif detset in ["70B", "070B"]:
        horns = [19, 21, 22]
        instrument = "LFI"
    elif isinstance(detset, str) and detset.upper() == "LFI":
        detectors.extend(list_planck(30, good=good))
        detectors.extend(list_planck(44, good=good))
        detectors.extend(list_planck(70, good=good))
        return detectors

    # --- HFI channel selections ------------------------------------------------
    elif detset in (100, "100", "100GHz"):
        psb_horns = range(1, 5)
        swb_horns: Sequence[int] = []
        if subset == 1:
            psb_horns = [1, 4]
        elif subset == 2:
            psb_horns = [2, 3]
        instrument = "HFI"
        freq = "100-"
    elif detset == "100A":
        psb_horns = [1, 4]
        swb_horns = []
        instrument = "HFI"
        freq = "100-"
    elif detset == "100B":
        psb_horns = [2, 3]
        swb_horns = []
        instrument = "HFI"
        freq = "100-"
    elif detset in (143, "143", "143GHz"):
        psb_horns = np.arange(1, 5)
        if good:
            swb_horns = range(5, 8)
        else:
            swb_horns = range(5, 9)
        if subset == 1:
            psb_horns, swb_horns = [1, 3], []
        elif subset == 2:
            psb_horns, swb_horns = [2, 4], []
        elif subset == 3:
            psb_horns, swb_horns = [], [5, 6, 7]
        instrument = "HFI"
        freq = "143-"
    elif detset == "143A":
        psb_horns = [1, 3]
        swb_horns = [5, 7]
        instrument = "HFI"
        freq = "143-"
    elif detset == "143B":
        psb_horns = [2, 4]
        swb_horns = [6]
        instrument = "HFI"
        freq = "143-"
    elif detset in (217, "217", "217GHz"):
        psb_horns = np.arange(5, 9)
        swb_horns = np.arange(1, 5)
        if subset == 1:
            psb_horns, swb_horns = [5, 7], []
        elif subset == 2:
            psb_horns, swb_horns = [6, 8], []
        elif subset == 3:
            psb_horns, swb_horns = [], [1, 2, 3, 4]
        instrument = "HFI"
        freq = "217-"
    elif detset == "217A":
        psb_horns = [5, 7]
        swb_horns = [1, 3]
        instrument = "HFI"
        freq = "217-"
    elif detset == "217B":
        psb_horns = [6, 8]
        swb_horns = [2, 4]
        instrument = "HFI"
        freq = "217-"
    elif detset in (353, "353", "353GHz"):
        psb_horns = np.arange(3, 7)
        swb_horns = [1, 2, 7, 8]
        if subset == 1:
            psb_horns, swb_horns = [3, 5], []
        elif subset == 2:
            psb_horns, swb_horns = [4, 6], []
        elif subset == 3:
            psb_horns, swb_horns = [], [1, 2, 7, 8]
        instrument = "HFI"
        freq = "353-"
    elif detset == "353A":
        psb_horns = [3, 5]
        swb_horns = [1, 7]
        instrument = "HFI"
        freq = "353-"
    elif detset == "353B":
        psb_horns = [4, 6]
        swb_horns = [2, 8]
        instrument = "HFI"
        freq = "353-"
    elif detset in (545, "545", "545GHz"):
        psb_horns = []
        if good and not extend_545:
            swb_horns = [1, 2, 4]
        else:
            swb_horns = np.arange(1, 5)
        instrument = "HFI"
        freq = "545-"
    elif detset == "545A":
        psb_horns = []
        swb_horns = [1]
        instrument = "HFI"
        freq = "545-"
    elif detset == "545B":
        psb_horns = []
        swb_horns = [2, 4]
        instrument = "HFI"
        freq = "545-"
    elif detset in (857, "857", "857GHz"):
        psb_horns = []
        if good and not extend_857:
            swb_horns = [1, 2, 3]
        else:
            swb_horns = np.arange(1, 5)
        instrument = "HFI"
        freq = "857-"
    elif detset == "857A":
        psb_horns = []
        swb_horns = [1, 3]
        instrument = "HFI"
        freq = "857-"
    elif detset == "857B":
        psb_horns = []
        swb_horns = [2, 4]
        instrument = "HFI"
        freq = "857-"
    elif isinstance(detset, str) and detset.upper() == "HFI":
        detectors.extend(list_planck(100, good=good, extend_857=extend_857))
        detectors.extend(list_planck(143, good=good, extend_857=extend_857))
        detectors.extend(list_planck(217, good=good, extend_857=extend_857))
        detectors.extend(list_planck(353, good=good, extend_857=extend_857))
        detectors.extend(list_planck(545, good=good, extend_857=extend_857))
        detectors.extend(list_planck(857, good=good, extend_857=extend_857))
        return detectors
    elif isinstance(detset, str) and detset.upper() == "PLANCK":
        detectors.extend(list_planck("LFI", good=good, extend_857=extend_857))
        detectors.extend(list_planck("HFI", good=good, extend_857=extend_857))
        return detectors
    elif isinstance(detset, str) and detset.upper() == "ROWS":
        # Row-grouping of detectors; kept as in the original code
        return [
            ["LFI27M", "LFI27S", "LFI28M", "LFI28S"],
            ["LFI24M", "LFI24S"],
            ["LFI25M", "LFI25S", "LFI26M", "LFI26S"],
            ["LFI18M", "LFI23S"],
            ["LFI19M", "LFI22S"],
            ["LFI20M", "LFI21S"],
            ["100-1a", "100-1b", "100-4a", "100-4b"],
            ["100-2a", "100-2b", "100-3a", "100-3b"],
            ["143-1a", "143-1b", "143-3a", "143-3b"],
            ["143-2a", "143-2b", "143-4a", "143-4b"],
            ["143-5", "143-7"],
            ["143-6"],
            ["217-1", "217-3"],
            ["217-2", "217-4"],
            ["217-5a", "217-5b", "217-7a", "217-7b"],
            ["217-6a", "217-6b", "217-8a", "217-8b"],
            ["353-1", "353-7"],
            ["353-3a", "353-3b", "353-5a", "353-5b"],
            ["353-4a", "353-4b", "353-6a", "353-6b"],
            ["353-2", "353-8"],
            ["545-1"],
            ["545-2", "545-4"],
            ["857-1", "857-3"],
            ["857-2"],
        ]
    else:
        # Single detectors and horns
        lfidets = list_planck("LFI")
        hfidets = list_planck("HFI")
        if detset in lfidets or detset in hfidets:
            return [detset]  # type: ignore[list-item]
        if isinstance(detset, str) and detset + "M" in lfidets:
            return [detset + "M", detset + "S"]  # type: ignore[return-value]
        if isinstance(detset, str) and detset + "a" in hfidets:
            return [detset + "a", detset + "b"]  # type: ignore[return-value]
        # All other cases
        print("ERROR: unknown detector set: ", detset)
        return -1

    # Build detector name list for LFI / HFI cases above
    if instrument == "LFI":
        for horn in horns:
            for arm in ["S", "M"]:
                detectors.append("LFI" + str(horn) + arm)
    elif instrument == "HFI":
        for horn in psb_horns:
            for arm in ["a", "b"]:
                detectors.append(freq + str(horn) + arm)
        for horn in swb_horns:
            detectors.append(freq + str(horn))

    return detectors


# ----------------------------------------------------------------------
# Detector weights
# ----------------------------------------------------------------------

#: Map-level detector weights (for example from white-noise levels).
detector_weights: Dict[str, float] = {
    # 30 GHz
    "LFI27": 0.40164e06,
    "LFI28": 0.36900e06,
    # 44 GHz
    "LFI24": 0.12372e06,
    "LFI25": 0.14049e06,
    "LFI26": 0.11233e06,
    # 70 GHz
    "LFI18": 53650.0,
    "LFI19": 42141.0,
    "LFI20": 36579.0,
    "LFI21": 50355.0,
    "LFI22": 49363.0,
    "LFI23": 47966.0,
    # 100 GHz
    "100-1": 0.76343e06,
    "100-2": 0.12661e07,
    "100-3": 0.10631e07,
    "100-4": 0.10532e07,
    # 143 GHz
    "143-1": 0.16407e07,
    "143-2": 0.18577e07,
    "143-3": 0.16439e07,
    "143-4": 0.14458e07,
    "143-5": 0.27630e07,
    "143-6": 0.26942e07,
    "143-7": 0.28599e07,
    # 217 GHz
    "217-1": 0.11058e07,
    "217-2": 0.10261e07,
    "217-3": 0.10958e07,
    "217-4": 0.10593e07,
    "217-5": 0.67318e06,
    "217-6": 0.71092e06,
    "217-7": 0.76576e06,
    "217-8": 0.71226e06,
    # 353 GHz
    "353-1": 0.12829e06,
    "353-2": 0.13475e06,
    "353-3": 48067.0,
    "353-4": 42187.0,
    "353-5": 56914.0,
    "353-6": 25293.0,
    "353-7": 87730.0,
    "353-8": 74453.0,
    # 545 GHz
    "545-1": 4475.5,
    "545-2": 5540.3,
    "545-4": 4321.0,
    # 857 GHz
    "857-1": 6.8895,
    "857-2": 6.3108,
    "857-3": 6.5964,
    "857-4": 3.6785,
}


# ----------------------------------------------------------------------
# QuickPol filename helper
# ----------------------------------------------------------------------


def qp_file(
    outdir: str,
    dets: Sequence[str],
    lmax: Optional[int] = 2000,
    smax: Optional[int] = 6,
    angle_shift: Optional[float] = 0.0,
    full: Optional[bool] = True,
    pconv: str = "cmbfast",
    force_det: Optional[str] = None,
    release: Optional[str] = None,  # kept for compatibility, not used
    rhobeam: Optional[str] = None,
    rhohit: Optional[str] = None,
) -> str:
    """
    Construct the standard QuickPol beam matrix filename.

    Parameters
    ----------
    outdir : str
        Directory where the file will live. Created if it does not exist.
    dets : sequence of str
        Two-element sequence with detector identifiers, e.g. ['143-1a', '143-1b'].
    lmax : int, optional
        Maximum multipole ℓ used in QuickPol (default: 2000).
    smax : int, optional
        Maximum spin moment used in the expansion (default: 6).
    angle_shift : float, optional
        Additional angle label (in degrees) for the filename. Used only in the
        name, not in any internal computation.
    full : bool, optional
        If True, label the file as containing full beam information.
    pconv : str, optional
        Power spectrum convention label (e.g. 'cmbfast').
    force_det : str, optional
        If not None, add a '_FD<force_det>' tag to the filename.
    release : str, optional
        Reserved for data-release dependent tags (currently unused).
    rhobeam : {'Ideal', 'IMO'}, optional
        Beam mismatch model label. 'Ideal' -> no tag, 'IMO' -> '_rbIMO'.
    rhohit : {'Ideal', 'IMO'}, optional
        Hit mismatch model label. 'Ideal' -> no tag, 'IMO' -> '_rhIMO'.

    Returns
    -------
    str
        Full path to the QuickPol beam matrix file (NPZ).
    """
    # Ensure output directory exists
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    lmax_def = 3000
    if lmax is None:
        lmax = lmax_def
    if smax is None:
        smax = 6
    if angle_shift is None:
        angle_shift = 0.0
    if full is None:
        full = True

    if force_det is not None:
        sfd = "_FD%s" % (force_det)
    else:
        sfd = ""

    if rhobeam == "Ideal":
        srb = ""
    elif rhobeam == "IMO":
        srb = "_rbIMO"
    elif rhobeam is None:
        srb = ""
    else:
        raise RuntimeError(f"Unknown rhobeam: {rhobeam}")

    if rhohit == "Ideal":
        srh = ""
    elif rhohit == "IMO":
        srh = "_rhIMO"
    elif rhohit is None:
        srh = ""
    else:
        raise RuntimeError(f"Unknown rhohit: {rhohit}")

    angst = "%+03d" % (angle_shift)
    angst = angst.replace("+180", "180")
    angst = angst.replace("-180", "180")
    angst = angst.replace("+00", "000")

    fz = os.path.join(
        outdir,
        "beam_matrix_{}x{}_l{}_s{}_A{}_{}_{}{}{}{}.npz".format(
            dets[0],
            dets[1],
            str(lmax),
            str(smax),
            angst,
            pconv,
            str(int(bool(full))),
            sfd,
            srb,
            srh,
        ),
    )

    return fz
