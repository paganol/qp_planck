"""
QuickPol beam window → FITS converter for Planck NPIPE.

DISCLAIMER
----------
This file contains code adapted from the Planck NPIPE pipeline:
    https://github.com/planck-npipe/toast-npipe/tree/master

The original implementation relied on TOAST and internal Planck tooling.
This version removes TOAST dependencies and is intended to be used as
a standalone step to convert QuickPol NPZ outputs into FITS window
function files (B_ell and W_ell).
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

"""
    main
      +--q2f
           +--clobber
           +--my_mwrfits
"""

import datetime
import os

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

import astropy.io.fits as pyfits
import healpy as hp

from utilities import list_planck, qp_file

outdir = '../quickpol_output'
indir = '../quickpol_output'
smax = 6
docross = True
blfile = True
wlfile = True
blTEBfile = True
overwrite = False
release = 'npipe6v20'
full = False  # False : Only sample a small fraction of the pixels
do_plot = False

NO_COLOR = '\x1b[0m'
GREEN_COLOR = '\x1b[32;11m'
RED_COLOR = '\x1b[31;01m'
BLUE_COLOR = '\x1b[34;11m'
BOLD = '\x1b[1;01m'

t1 = np.array([  # sym
    ['TT', 'TE', 'TB'],
    ['TE', 'EE', 'EB'],
    ['TB', 'EB', 'BB']])
t2 = np.array([  # non-sym
    ['TT', 'TE', 'TB'],
    ['ET', 'EE', 'EB'],
    ['BT', 'BE', 'BB']])
kk = [[0, 0], [1, 1], [2, 2], [0, 1], [0, 2], [1, 2], [1, 0], [2, 0], [2, 1]]
t3 = [t2[k1, k2] for k1, k2 in kk]


# ==============================================================================


def fit_gauss(bl):
    """
    Fit a Gaussian beam to the provided beam window function.

    Parameters
    ----------
    bl : array_like
        Beam window function B_ell (1D array over ell).

    Returns
    -------
    gaussbeam : ndarray
        Best-fit Gaussian beam evaluated at all ell.
    sigma : float
        Best-fit Gaussian width in radians.
    """
    ell = np.arange(bl.size)

    def gaussbeam(ell, sigma):
        return np.exp(-0.5 * ell * (ell + 1) * sigma ** 2)

    def resid(p, ell, bl):
        sigma = p[0]
        return bl - gaussbeam(ell, sigma)

    p0 = [np.radians(0.5)]
    result = scipy.optimize.least_squares(
        resid, p0, method='lm', args=(ell, bl), max_nfev=10000
    )
    if not result.success:
        raise RuntimeError(f'Gaussian fitting failed: {result.message}')
    sigma = result.x[0]
    return gaussbeam(ell, sigma), sigma


# ------------------------------------------------------------------------------


def clobber(filename, overwrite):
    """
    Decide whether to write a file, respecting an overwrite flag.

    Parameters
    ----------
    filename : str
        Path of the file to write.
    overwrite : bool
        Whether to overwrite if the file already exists.

    Returns
    -------
    bool
        True if the caller should write the file, False otherwise.
    """
    write = True
    if os.path.exists(filename):
        if overwrite:
            print(prefix, f'{RED_COLOR}Overwriting {filename}{NO_COLOR}',
                  flush=True)
        else:
            print(prefix,
                  f'{BLUE_COLOR}{filename} already exists. Skip{NO_COLOR}',
                  flush=True)
            write = False
    return write


# ------------------------------------------------------------------------------


def my_mwrfits(
    filename,
    data,
    colnames=None,
    keys=None,
    bintable=False,
    ftype=None,
    extnames=None,
    origin=None,
    dets=None,
):
    """Write columns to a FITS file in one or more table extensions.

    Parameters
    ----------
    filename : str
        The FITS file name.
    data : list[list[np.ndarray]]
        A list (per extension) of lists (per column) of 1D arrays.
    colnames : list[list[str]]
        Column names per extension.
    keys : dict-like, optional
        Extra header keywords to write in each table header.
    bintable : bool, optional
        If True, write binary table HDUs; otherwise ASCII table HDUs.
    ftype : {'B', 'B_TEB', 'W'}, optional
        Type of window function (B_ell, B_ell T/E/B, or W_ell).
    extnames : list[str], optional
        Names of the FITS extensions.
    origin : list[str], optional
        Comment lines describing provenance (e.g. code, NPZ input).
    dets : list[str], optional
        Detector or detector-set identifiers.
    """
    hline = '----------------------------------------------------------------'
    if ftype == 'B':
        comments = [
            'Beam Window Function B(l)',
            'Compatible with Healpix (synfast, smoothing, ...) and PolSpice',
            'To be squared before applying to power spectrum',
            '  C_map(l) = C_sky(l) * B(l)^2 ',
        ]
    elif ftype == 'B_TEB':
        comments = [
            'Beam Window Functions B(l), for T, E and B',
            'Compatible with Healpix (synfast, smoothing, ...) and PolSpice',
            'To be squared before applying to power spectrum',
            '  C_TT_map(l) = C_TT_sky(l) * B_T(l)^2 ',
            '  C_EE_map(l) = C_EE_sky(l) * B_E(l)^2 ',
            '  C_BB_map(l) = C_BB_sky(l) * B_B(l)^2 ',
        ]
    elif ftype == 'W':
        comments = [
            'Beam Window Functions W(l) = B(l)^2',
            'Applies directly to power spectrum',
            '  C_map(l) = C_sky(l) * W(l) ',
            'Includes cross-talk terms',
        ]
    else:
        comments = []

    # ---- primary header -----
    hdu = pyfits.PrimaryHDU(None)
    hhu = hdu.header.set
    hhc = hdu.header.add_comment
    fdate = datetime.datetime.now().strftime('%Y-%m-%d')
    hhu('DATE', fdate, comment='Creation date (CCYY-MM-DD) of FITS header')
    if extnames is not None:
        nx = len(extnames)
        hhu('NUMEXT', nx, 'Number of extensions')
        for xt in range(nx):
            hhu(f'XTNAME{xt+1}', extnames[xt],
                f'Name of extension #{xt+1}')

    hhc(hline)
    for mycom in comments:
        hhc(mycom)
    if origin is not None:
        for myor in origin:
            hhc(myor)
        hhc(hline)
    if dets is not None:
        for id, det in enumerate(dets):
            hhu(f'DET{id+1}', det, 'Detector (set)')

    hdulist = pyfits.HDUList([hdu])

    # ---- other HDUs : tables ----
    getformat = hp.fitsfunc.getformat

    for xt in range(len(data)):
        cols = []
        for line in range(len(data[xt])):
            namei = colnames[xt][line]
            array = data[xt][line]
            if bintable:
                nt = len(array)  # total length
                repeat = nt  # length / cell
                fmt = str(repeat) + getformat(array)
                array = np.reshape(array, (nt // repeat, repeat))
            else:
                fmt = getformat(array)
            cols.append(pyfits.Column(name=namei,
                                      format=fmt,
                                      array=array))
            if bintable:
                tbhdu = pyfits.BinTableHDU.from_columns(cols)
            else:
                tbhdu = pyfits.TableHDU.from_columns(cols)

        if extnames is not None:
            tbhdu.name = extnames[xt]

        ncols = len(cols)
        tbhdu.header['MAX-LPOL'] = (len(data[xt][0]) - 1,
                                    'Maximum L multipole')
        tbhdu.header['POLAR'] = (ncols > 1)
        tbhdu.header['BCROSS'] = (ncols > 4)
        tbhdu.header['ASYMCL'] = (ncols > 6)

        tbhdu.header.add_comment(hline)
        for mycom in comments:
            tbhdu.header.add_comment(mycom)
        if origin is not None:
            for myor in origin:
                tbhdu.header.add_comment(myor)
        tbhdu.header.add_comment(hline)

        if isinstance(keys, dict):
            for k, v in list(keys.items()):
                tbhdu.header[k] = (v)

        hdulist.append(tbhdu)

    # write the file
    hdulist.writeto(filename, overwrite=True)

    # checking out the file
    try:
        _ = pyfits.getdata(filename)
        _ = hp.mrdfits(filename)
        print(prefix, f'{GREEN_COLOR} checking out {filename}{NO_COLOR}',
              flush=True)
    except Exception as err:
        raise RuntimeError(f'Failed to load {filename}: {err}') from err


# -----------------------------------------------------------------------------


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
    if detset.startswith('0') or detset.startswith('LFI'):
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
    if '545' in detset or '857' in detset or 'LFI' in detset or '-' in detset:
        pol = False
    else:
        pol = True
    return pol


def q2f(indir, outdir, dets, smax, release=None, full=True, blfile=True,
        blTEBfile=True, wlfile=True, overwrite=True, do_plot=False):
    """
    Convert QuickPol NPZ output into FITS B_ell and W_ell window function files.

    Parameters
    ----------
    indir : str
        Directory containing the QuickPol NPZ files.
    outdir : str
        Directory where the FITS window function files will be written.
    dets : (str, str)
        Detector or detector-set pair (e.g. ('143A', '143A')).
    smax : int
        Maximum spin used in the QuickPol computation (just passed through).
    release : str, optional
        Release tag used in output filenames, e.g. 'npipe6v20'.
    full : bool, optional
        Whether the QuickPol output corresponds to full sampling (used only in
        the NPZ filename construction).
    blfile : bool, optional
        Write scalar B_ell (temperature) FITS file.
    blTEBfile : bool, optional
        Write B_ell for T, E, B if polarization is available.
    wlfile : bool, optional
        Write W_ell matrices in a binary table.
    overwrite : bool, optional
        Overwrite existing files.
    do_plot : bool, optional
        If True, save simple PNG plots of the window functions.
    """
    pconv = 'cmbfast'
    angle_shift = 0
    force_det = None
    rhobeam = 'IMO'
    rhohit = 'IMO'

    lmax = min(detset2lmax(dets[0]), detset2lmax(dets[1]))
    pol = detset2pol(dets[0]) and detset2pol(dets[1])

    fz = qp_file(indir, dets, lmax=lmax, smax=smax, angle_shift=angle_shift,
                 full=full, pconv=pconv, force_det=force_det, release=release,
                 rhobeam=rhobeam, rhohit=rhohit)
    print(prefix, '--------------------')
    print(prefix, fz, flush=True)
    try:
        dz1 = np.load(fz, allow_pickle=True)
    except Exception:
        print(prefix, f'{fz} not found', flush=True)
        return

    f32 = np.float32
    bm1 = dz1['beam_mat'].tolist()
    TT = f32(bm1['TT'])
    renorm = TT[0, 0, 0]
    TT /= renorm
    EE = f32(bm1['EE']) / renorm
    BB = f32(bm1['BB']) / renorm
    TE = f32(bm1['TE']) / renorm
    print(prefix,
          f'{BLUE_COLOR} Renorm-1 = {renorm - 1}{NO_COLOR}',
          flush=True)
    wtt = TT[0:lmax + 1, 0, 0]
    bl = np.sqrt(np.abs(wtt)) * np.sign(wtt)
    imin = np.argmin(bl)
    imax = np.argmax(bl)
    wee = EE[0:lmax + 1, 1, 1]
    wbb = BB[0:lmax + 1, 2, 2]
    bl_E = np.sqrt(np.abs(wee)) * np.sign(wee)
    bl_B = np.sqrt(np.abs(wbb)) * np.sign(wbb)

    ineg = np.where(bl < 0)[0]
    print(prefix, 'Max = ', bl[imax], imax)
    if len(ineg) > 0:
        print(prefix,
              f'{RED_COLOR} Neg = {ineg[0]} {ineg[-1]}{NO_COLOR}',
              flush=True)
    print(prefix, 'Min = ', bl[imin], imin, flush=True)

    fitsfile_T = os.path.join(
        outdir, f'Bl_{release}_{dets[0]}x{dets[1]}.fits')
    fitsfile_TEB = os.path.join(
        outdir, f'Bl_TEB_{release}_{dets[0]}x{dets[1]}.fits')
    fitsfile_W = os.path.join(
        outdir, f'Wl_{release}_{dets[0]}x{dets[1]}.fits')

    fdate = datetime.datetime.now().strftime('%Y-%m-%d')
    origin = ['Adapted from', fz, f'by {__file__} on {fdate}']

    gaussbeam, sigma = fit_gauss(bl)
    fwhm = np.abs(np.degrees(sigma) * 60 * np.sqrt(8.0 * np.log(2.0)))

    # T B(l)
    if blfile and clobber(fitsfile_T, overwrite):
        extnames = ['WINDOW FUNCTION']
        my_mwrfits(
            fitsfile_T,
            [[bl]],
            colnames=[['TEMPERATURE']],
            bintable=False,
            ftype='B',
            extnames=extnames,
            origin=origin,
            dets=dets,
        )
        if do_plot:
            hdulist = pyfits.open(fitsfile_T)
            plt.figure()
            plt.gca().set_title(f'{release} {dets[0]} x {dets[1]}')
            plt.semilogy(hdulist[1].data.field(0), label='T')
            ylim = [1e-8, 2]
            plt.plot(gaussbeam, label=f"{fwhm:.2f}' FWHM")
            plt.gca().set_ylim(ylim)
            plt.legend(loc='best')
            fn_plot = fitsfile_T.replace('.fits', '.png')
            plt.savefig(fn_plot)
            print(prefix, 'Plot saved in', fn_plot, flush=True)
            plt.close()
            hdulist.close()

    # T, E, B B(l)
    if blTEBfile and clobber(fitsfile_TEB, overwrite) and pol:
        extnames = ['WINDOW FUNCTIONS']
        my_mwrfits(
            fitsfile_TEB,
            [[bl, bl_E, bl_B]],
            colnames=[['T', 'E', 'B']],
            bintable=False,
            ftype='B_TEB',
            extnames=extnames,
            origin=origin,
            dets=dets,
        )
        if do_plot:
            hdulist = pyfits.open(fitsfile_TEB)
            plt.figure()
            plt.gca().set_title(f'{release} {dets[0]} x {dets[1]}')
            for i, lab in enumerate('TEB'):
                plt.semilogy(hdulist[1].data.field(i), label=lab)
            ylim = [1e-8, 2]
            plt.plot(gaussbeam, label=f"{fwhm:.2f}' FWHM")
            plt.gca().set_ylim(ylim)
            plt.legend(loc='best')
            fn_plot = fitsfile_TEB.replace('.fits', '.png')
            plt.savefig(fn_plot)
            print(prefix, 'Plot saved in', fn_plot, flush=True)
            plt.close()
            hdulist.close()

    # W(l)
    if wlfile and clobber(fitsfile_W, overwrite) and pol:
        extnames = ['TT', 'EE', 'BB', 'TE']
        data = [
            [TT[0:, k1, k2] for k1, k2 in kk],
            [EE[0:, k1, k2] for k1, k2 in kk],
            [BB[0:, k1, k2] for k1, k2 in kk],
            [TE[0:, k1, k2] for k1, k2 in kk],
        ]
        colnames = [
            [extnames[0] + '_2_' + c for c in t3],
            [extnames[1] + '_2_' + c for c in t3],
            [extnames[2] + '_2_' + c for c in t3],
            [extnames[3] + '_2_' + c for c in t3],
        ]
        my_mwrfits(
            fitsfile_W,
            data,
            colnames=colnames,
            bintable=True,
            ftype='W',
            extnames=extnames,
            origin=origin,
            dets=dets,
        )
        if do_plot:
            hdulist = pyfits.open(fitsfile_W)
            plt.figure(figsize=[18, 12])
            plt.gca().set_title(f'{release} {dets[0]} x {dets[1]}')
            for ifield, field in enumerate(['TT_2_TE', 'TT_2_EE', 'TT_2_BB']):
                plt.subplot(2, 2, 1 + ifield)
                plt.plot(hdulist[1].data.field(field).flatten(), label=field)
                plt.legend(loc='best')
                plt.gca().axhline(0, color='k')
            for ifield, field in enumerate(['EE_2_BB']):
                plt.subplot(2, 2, 4 + ifield)
                plt.plot(hdulist[2].data.field(field).flatten(), label=field)
                plt.legend(loc='best')
                plt.gca().axhline(0, color='k')
            fn_plot = fitsfile_W.replace('.fits', '.png')
            plt.savefig(fn_plot)
            print(prefix, 'Plot saved in', fn_plot, flush=True)
            plt.close()
            hdulist.close()


# ------------------------------------------------------------------------------
# -------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if __name__ == '__main__':

    freqs = [30, 44, 70, 100, 143, 217, 353, 545, 857]

    detsets = []
    for suffix in ['GHz', 'A', 'B']:
        for freq in freqs:
            detset = '{:03}{}'.format(freq, suffix)
            detsets.append(detset)

    detsetpairs = []

    # Full frequency and detector set auto and cross spectra
    for idetset1, detset1 in enumerate(detsets):
        for idetset2, detset2 in enumerate(detsets):
            # if idetset2 < idetset1:
            #     continue
            # No cross spectra between full frequency and detsets
            if detset1.endswith('GHz') and detset2[-1] in 'AB':
                continue
            if detset2.endswith('GHz') and detset1[-1] in 'AB':
                continue
            detsetpairs.append((detset1, detset2))

    # Single detector and single horn auto spectra
    for det in list_planck('Planck'):
        # Single detector
        detsetpairs.append((det, det))
        if det[-1] in 'aM':
            # Single horn
            horn = det[:-1]
            detsetpairs.append((horn, horn))

    for ipair, detsetpair in enumerate(detsetpairs):
        if ipair % ntask != rank:
            continue
        q2f(outdir, indir, detsetpair, smax, release=release, full=full,
            blfile=blfile, blTEBfile=blTEBfile, wlfile=wlfile,
            overwrite=overwrite, do_plot=do_plot)

