# -*- coding: utf-8 -*-
"""
Created on Apr 11 2019

@author: fmillour
"""
from shutil import copyfile
from pathlib import Path
from typing import Optional

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


def read_fits(fits_file: Path,
              fits_file_replace: Optional[Path] = None,
              rvis_type: Optional[bool] = None) -> fits.HDUList:
    """Reads a FITS file.

    Parameters
    ----------
    fits_file : pathlib.Path
        The file to be read.
    fits_file_replace : pathlib.Path
        A file to replace the orginal one if it is None.
        Will create arrays with only zeros.

    Returns
    -------
    fits.HDUList
        The HDUList of the FITS file.
    """
    if fits_file in [None, "", " "]:
        with fits.open(fits_file_replace, "readonly") as hdul:
            t3phi = np.zeros_like(hdul["OI_T3"].data["T3PHI"])
            visphi = np.zeros_like(hdul["OI_VIS"].data["VISPHI"])
            visamp = np.zeros_like(hdul["OI_VIS"].data["VISAMP"])
            vis2data = np.zeros_like(hdul["OI_VIS2"].data["VIS2DATA"])
            vis_type = np.zeros_like(hdul["OI_VIS"].header["AMPTYP"])
    else:
        with fits.open(fits_file, "readonly") as hdul:
            t3phi = hdul["OI_T3"].data["T3PHI"]*u.deg.to(u.rad)
            visphi = hdul["OI_VIS"].data["VISPHI"]*u.deg.to(u.rad)
            visamp = hdul["OI_VIS"].data["VISAMP"]
            vis2data = hdul["OI_VIS2"].data["VIS2DATA"]
            wavelength = hdul["OI_WAVELENGTH"].data["EFF_WAVE"]
            vis_type = hdul["OI_VIS"].header["AMPTYP"]
    if rvis_type:
        return t3phi, visphi, visamp, vis2data, wavelength, vis_type
    return t3phi, visphi, visamp, vis2data


def compact_exposures(nrepeat: int,
                      vis: np.ndarray, vis2: np.ndarray,
                      sin_dp: np.ndarray, cos_dp: np.ndarray,
                      sin_t3phi: np.ndarray, cos_t3phi: np.ndarray):
    """Stores multiple exposures in the first six rows."""
    if nrepeat > 1:
        for i in np.arange(nrepeat - 1):
            for j in np.arange(6):
                vis[j, :] += vis[(i + 1) * 6 + j, :]
                vis2[j, :] += vis2[(i + 1) * 6 + j, :]
                sin_dp[j, :] += sin_dp[(i + 1) * 6 + j, :]
                cos_dp[j, :] += cos_dp[(i + 1) * 6 + j, :]

            for j in np.arange(4):
                sin_t3phi[j, :] += sin_t3phi[(i + 1) * 4 + j, :]
                cos_t3phi[j, :] += cos_t3phi[(i + 1) * 4 + j, :]
    return vis, vis2, sin_dp, cos_dp, sin_t3phi, cos_t3phi


def calib_BCD(iifile: Path, iofile: Path,
              oifile: Path, oofile: Path,
              outputfile: Optional[Path] = None,
              lim: Optional[int] = 180,
              plot: Optional[bool] = True) -> None:
    """Calibrates the different exposures of the Beam Commuting
    Device (BCD).

    Parameters
    ----------
    iifile : pathlib.Path
        The in-in file.
    iofile : pathlib.Path
        The in-out file.
    oifile : pathlib.Path
        The out-in file.
    oofile : pathlib.Path
        The out-out file.
    outputfile : pathlib.Path, optional
        The name of the output file, by default Path.cwd() / "toto.fits"
    lim : int, optional
        The number of degrees to limit the angle, by default 180
    plot : bool, optional
        Whether to plot the results, by default True
    """
    if outputfile is None:
        outputfile = Path.cwd() / "toto.fits"
    copyfile(iifile, outputfile)

    iifile, oofile = map(Path, [iifile, oofile])
    iit3p, iidp, iiva, iiv2, iitwl, vis_type = read_fits(iifile,
                                                         rvis_type=True)
    oot3p, oodp, oova, oov2 = read_fits(oofile)
    iot3p, iodp, iova, iov2 = read_fits(iofile, iifile)
    oit3p, oidp, oiva, oiv2 = read_fits(oifile, iifile)

    sin_iidp, sin_iit3p = np.sin(iidp), np.sin(iit3p)
    sin_iodp, sin_iot3p = np.sin(iodp), np.sin(iot3p)
    sin_oidp, sin_oit3p = np.sin(oidp), np.sin(oit3p)
    sin_oodp, sin_oot3p = np.sin(oodp), np.sin(oot3p)
    cos_iidp, cos_iit3p = np.cos(iidp), np.cos(iit3p)
    cos_iodp, cos_iot3p = np.cos(iodp), np.cos(iot3p)
    cos_oidp, cos_oit3p = np.cos(oidp), np.cos(oit3p)
    cos_oodp, cos_oot3p = np.cos(oodp), np.cos(oot3p)

    # NOTE: Addup the different exposures
    nwlen = np.shape(iit3p)[1]
    nrepeatii = int(np.shape(iiv2)[0] / 6)
    nrepeatio = int(np.shape(iov2)[0] / 6)\
        if iofile not in [None, " ", ""] else 0
    nrepeatoi = int(np.shape(oiv2)[0] / 6)\
        if oifile not in [None, " ", ""] else 0
    nrepeatoo = int(np.shape(oov2)[0] / 6)

    # NOTE: Store multiple exposures data into the first 6 rows
    iiva, iiv2, sin_iidp, cos_iidp, sin_iit3p, cos_iit3p = compact_exposures(
            nrepeatii, iiva, iiv2, sin_iidp, cos_iidp, sin_iit3p, cos_iit3p)

    oova, oov2, sin_oodp, cos_oodp, sin_oot3p, cos_oot3p = compact_exposures(
            nrepeatoo, oova, oov2, sin_oodp, cos_oodp, sin_oot3p, cos_oot3p)

    iova, iov2, sin_iodp, cos_iodp, sin_iot3p, cos_iot3p = compact_exposures(
            nrepeatio, iova, iov2, sin_iodp, cos_iodp, sin_iot3p, cos_iot3p)

    oiva, oiv2, sin_oidp, cos_oidp, sin_oit3p, cos_oit3p = compact_exposures(
            nrepeatoi, oiva, oiv2, sin_oidp, cos_oidp, sin_oit3p, cos_oit3p)

    # NOTE: Treat closure phases
    idx = np.array([[0, 0, 3, 3], [1, 2, 1, 2],
                    [2, 1, 2, 1], [3, 3, 0, 0]])
    sign = np.array([[1, -1, 1, -1], [1, 1, -1, -1],
                     [1, 1, -1, -1], [1, -1, 1, -1]])

    # NOTE: Initialize closure phase with same shape as input
    sin_avg = np.zeros((4, nwlen))
    cos_avg = np.zeros((4, nwlen))
    closfinal = np.zeros((4, nwlen))

    if plot:
        plt.figure(51)

    for i in np.arange(4):
        sin_avg[i, :] = (
            sign[i, 0] * sin_iit3p[idx[i, 0], :]
            + sign[i, 1] * sin_oit3p[idx[i, 1], :]
            + sign[i, 2] * sin_iot3p[idx[i, 2], :]
            + sign[i, 3] * sin_oot3p[idx[i, 3], :]
        ) / (nrepeatii + nrepeatoi + nrepeatio + nrepeatoo)

        cos_avg[i, :] = (
            1.0 * cos_iit3p[idx[i, 0], :]
            + 1.0 * cos_oit3p[idx[i, 1], :]
            + 1.0 * cos_iot3p[idx[i, 2], :]
            + 1.0 * cos_oot3p[idx[i, 3], :]
        ) / (nrepeatii + nrepeatoi + nrepeatio + nrepeatoo)

        closfinal[i, :] = np.arctan2(
                sin_avg[i, :], cos_avg[i, :])*u.rad.to(u.deg)

        if plot:
            plt.subplot(2, 2, i + 1)
            plt.plot(iitwl * 1e6, iit3p[i, :], label="II")
            plt.plot(iitwl * 1e6, oot3p[i, :], label="OO")
            plt.plot(iitwl * 1e6, iot3p[i, :], label="IO")
            plt.plot(iitwl * 1e6, oit3p[i, :], label="OI")
            plt.plot(iitwl * 1e6, closfinal[i, :])
            plt.ylabel("Closure phase")
            plt.legend()
            plt.ylim(-lim, lim)

    # NOTE: Treat differential phases
    idx = np.array(
        [
            [0, 0, 0, 0],
            [1, 1, 1, 1],
            [2, 3, 4, 5],
            [3, 2, 5, 4],
            [4, 5, 2, 3],
            [5, 4, 3, 2],
        ]
    )
    sign = np.array(
        [
            [1, -1, 1, -1],
            [1, 1, -1, -1],
            [1, 1, 1, 1],
            [1, 1, 1, 1],
            [1, 1, 1, 1],
            [1, 1, 1, 1],
        ]
    )

    sin_avg = np.zeros((6, nwlen))
    cos_avg = np.zeros((6, nwlen))
    dpfinal = np.zeros((6, nwlen))

    if plot:
        plt.figure(52)

    for i in np.arange(6):
        sin_avg[i, :] = (
            sign[i, 0] * sin_iidp[idx[i, 0], :]
            + sign[i, 1] * sin_oidp[idx[i, 1], :]
            + sign[i, 2] * sin_iodp[idx[i, 2], :]
            + sign[i, 3] * sin_oodp[idx[i, 3], :]
        ) / (nrepeatii + nrepeatoi + nrepeatio + nrepeatoo)
        cos_avg[i, :] = (
            1.0 * cos_iidp[idx[i, 0], :]
            + 1.0 * cos_oidp[idx[i, 1], :]
            + 1.0 * cos_iodp[idx[i, 2], :]
            + 1.0 * cos_oodp[idx[i, 3], :]
        ) / (nrepeatii + nrepeatoi + nrepeatio + nrepeatoo)
        dpfinal[i, :] = np.arctan2(sin_avg[i, :], cos_avg[i, :]) * 180.0 / np.pi
        if plot:
            plt.subplot(3, 2, i + 1)
            plt.plot(iitwl * 1e6, iidp[i, :], label="II")
            plt.plot(iitwl * 1e6, oodp[i, :], label="OO")
            plt.plot(iitwl * 1e6, iodp[i, :], label="IO")
            plt.plot(iitwl * 1e6, oidp[i, :], label="OI")
            plt.plot(iitwl * 1e6, dpfinal[i, :])
            plt.ylim(-lim, lim)
            plt.ylabel("Differential phase")
            plt.legend()

    # NOTE: Treat visamp (or correlated fluxes)
    vafinal = np.zeros((6, nwlen))

    if plot:
        plt.figure(53)
    for i in np.arange(6):
        vafinal[i, :] = (
            iiva[i, :] + oiva[idx[i, 1], :] + iova[idx[i, 2], :] + oova[idx[i, 3], :]
        ) / (nrepeatii + nrepeatoi + nrepeatio + nrepeatoo)
        if plot:
            plt.subplot(3, 2, i + 1)
            plt.plot(iitwl * 1e6, iiva[i, :], label="II")
            plt.plot(iitwl * 1e6, oova[i, :], label="OO")
            plt.plot(iitwl * 1e6, iova[i, :], label="IO")
            plt.plot(iitwl * 1e6, oiva[i, :], label="OI")
            plt.plot(iitwl * 1e6, vafinal[i, :])
            plt.ylabel("VISAMP")
            plt.legend()

    # NOTE: Treat visibilities
    v2final = np.zeros((6, nwlen))

    if plot:
        plt.figure(50)

    for i in np.arange(6):
        v2final[i, :] = (
            iiv2[i, :] + oiv2[idx[i, 1], :] + iov2[idx[i, 2], :] + oov2[idx[i, 3], :]
        ) / (nrepeatii + nrepeatoi + nrepeatio + nrepeatoo)
        if plot:
            plt.subplot(3, 2, i + 1)
            plt.plot(iitwl * 1e6, iiv2[i, :], label="II")
            plt.plot(iitwl * 1e6, oov2[i, :], label="OO")
            plt.plot(iitwl * 1e6, iov2[i, :], label="IO")
            plt.plot(iitwl * 1e6, oiv2[i, :], label="OI")
            plt.plot(iitwl * 1e6, v2final[i, :])
            plt.ylim(0, 1)
            plt.ylabel("Squared visibility")
            plt.legend()
    breakpoint()

    with fits.open(outputfile, mode="update") as hdul:
        hdul["OI_T3"].data = hdul["OI_T3"].data[0:4]
        hdul["OI_T3"].data["T3PHI"] = closfinal[0:4]
        hdul["OI_VIS"].data = hdul["OI_VIS"].data[0:6]
        hdul["OI_VIS"].data["VISPHI"] = dpfinal[0:6]
        hdul["OI_VIS"].data["VISAMP"] = vafinal
        hdul["OI_VIS2"].data = hdul["OI_VIS2"].data[0:6]
        hdul["OI_VIS2"].data["VIS2DATA"] = v2final

        if "correlated" == vis_type:
            hdul["OI_VIS"].header["AMPTYP"] = "correlated flux"

        del hdul[0].header["HIERARCH ESO INS BCD1 ID"]
        del hdul[0].header["HIERARCH ESO INS BCD2 ID"]
        del hdul[0].header["HIERARCH ESO INS BCD1 NAME"]
        del hdul[0].header["HIERARCH ESO INS BCD2 NAME"]
        hdul[0].header["HIERARCH ESO INS BCD1 ID"] = " "
        hdul[0].header["HIERARCH ESO INS BCD2 ID"] = " "
        hdul[0].header["HIERARCH ESO INS BCD1 NAME"] = " "
        hdul[0].header["HIERARCH ESO INS BCD2 NAME"] = " "
        hdul.flush()

    if plot:
        pass
        plt.show()
