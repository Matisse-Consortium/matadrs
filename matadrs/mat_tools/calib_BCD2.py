# -*- coding: utf-8 -*-
"""
Created on Apr 11 2019

@author: fmillour
"""
from pathlib import Path
from typing import List, Optional, Union
from shutil import copyfile

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


T3P_INDEX = np.array([[0, 0, 3, 3], [1, 2, 1, 2],
                      [2, 1, 2, 1], [3, 3, 0, 0]])
T3P_SIGN = np.array([[1, -1, 1, -1], [1, 1, -1, -1],
                     [1, 1, -1, -1], [1, -1, 1, -1]])


DP_INDEX = np.array(
        [
            [0, 0, 0, 0],
            [1, 1, 1, 1],
            [2, 3, 4, 5],
            [3, 2, 5, 4],
            [4, 5, 2, 3],
            [5, 4, 3, 2],
            ]
        )

DP_SIGN = np.array(
        [
            [1, -1, 1, -1],
            [1, 1, -1, -1],
            [1, 1, 1, 1],
            [1, 1, 1, 1],
            [1, 1, 1, 1],
            [1, 1, 1, 1],
            ]
        )


def read_fits(fits_file: Path) -> fits.HDUList:
    """Reads a FITS file.

    Parameters
    ----------
    fits_file : pathlib.Path
        The file to be read.

    Returns
    -------
    fits.HDUList
        The HDUList of the FITS file.
    """
    with fits.open(fits_file, "readonly") as hdul:
        t3phi = hdul["OI_T3"].data["T3PHI"]*u.deg.to(u.rad)
        visphi = hdul["OI_VIS"].data["VISPHI"]*u.deg.to(u.rad)
        visamp = hdul["OI_VIS"].data["VISAMP"]
        vis2data = hdul["OI_VIS2"].data["VIS2DATA"]
        wavelength = hdul["OI_WAVELENGTH"].data["EFF_WAVE"]
        vis_type = hdul["OI_VIS"].header["AMPTYP"]
    return t3phi, visphi, visamp, vis2data, wavelength, vis_type


def compact_exposures(array: np.ndarray, compact_factor: int):
    """Stores multiple exposures in the first n rows given by the
    compact factor.

    Parameters
    ----------
    array : numpy.ndarray
        The array to be compacted.
    compact_factor : int
        The factor by which to compact the array.

    Returns
    -------
    compacted_array : numpy.ndarray
    """
    compact_array = array.copy()
    for i in range(array[0].shape[0]//compact_factor):
        for j in np.arange(compact_factor):
            compact_array[j] += array[(i + 1) * compact_factor + j]
    return compact_array


def plot(index: int, array: List[np.ndarray],
         final_array: np.ndarray,
         wavelength: np.ndarray, label: str,
         lim: Union[int, List[int]]) -> None:
    """Plots the array."""
    labels = ["II", "OO", "IO", "OI"]\
        if len(array) == 4 else ["II", "OO"]
    wavelength *= 1e6
    plt.subplot(len(array)//2, 2)
    for index, phase in enumerate(array):
        plt.plot(wavelength, phase, label=labels[index])
        plt.plot(wavelength, final_array[index])
    plt.ylabel(label)
    plt.ylim([-lim, lim] if isinstance(lim, int) else lim)
    plt.legend()


def calibrate_closure_phases(
        t3p: List[np.ndarray], wavelength: np.ndarray,
        file_factor: int, lim: Optional[int] = 180,
        do_plot: Optional[bool] = True) -> np.ndarray:
    """Calibrates the closure phases."""
    t3p = [compact_exposures(elem, compact_factor=4) for elem in t3p]
    sin_t3p, cos_t3p = zip(*[(np.sin(elem), np.cos(elem)) for elem in t3p])

    sin_avg = np.zeros((4, t3p[0].shape[1]))
    cos_avg = np.zeros((4, t3p[0].shape[1]))
    closfinal = np.zeros((4, t3p[0].shape[1]))

    if do_plot:
        plt.figure(51)

    for index in range(4):
        if len(t3p) == 4:
            sin_avg[index] = sum([T3P_SIGN[index, i]*elem[T3P_INDEX[index, i]]
                                  for i, elem in enumerate(sin_t3p)])
            cos_avg[index] = sum([elem[T3P_INDEX[index, i]]
                                  for i, elem in enumerate(cos_t3p)])
        else:
            indices = [0, 3]
            sin_avg[index] = sum([T3P_SIGN[index, i]*elem[T3P_INDEX[index, i]]
                                  for i, elem in zip(indices, sin_t3p)])
            cos_avg[index] = sum([elem[T3P_INDEX[index, i]]
                                  for i, elem in zip(indices, cos_t3p)])

        closfinal[index] = np.arctan2(
                sin_avg[index], cos_avg[index])*u.rad.to(u.deg)/file_factor

        if do_plot:
            plot(index, t3p, closfinal, wavelength,
                 "Closure phase", lim, do_plot)
    return closfinal


def calibrate_differential_phase(
        dp: List[np.ndarray], wavelength: np.ndarray,
        file_factor: int, lim: Optional[int] = 180,
        do_plot: Optional[bool] = True) -> np.ndarray:
    """Calibrates the differential phase."""
    dp = [compact_exposures(elem, compact_factor=4) for elem in dp]
    sin_dp, cos_dp = zip(*[(np.sin(elem), np.cos(elem)) for elem in dp])

    sin_avg = np.zeros((6, dp[0].shape[1]))
    cos_avg = np.zeros((6, dp[0].shape[1]))
    dpfinal = np.zeros((6, dp[0].shape[1]))

    if do_plot:
        plt.figure(52)

    for index in range(6):
        if len(dp) == 4:
            sin_avg[index] = sum([DP_SIGN[index, i]*elem[DP_INDEX[index, i]]
                                  for i, elem in enumerate(sin_dp)])
            cos_avg[index] = sum([elem[DP_INDEX[index, i]]
                                 for i, elem in enumerate(cos_dp)])
        else:
            indices = [0, 3]
            sin_avg[index] = sum([DP_SIGN[index, i]*elem[DP_INDEX[index, i]]
                                  for i, elem in zip(indices, sin_dp)])
            cos_avg[index] = sum([elem[DP_INDEX[index, i]]
                                 for i, elem in zip(indices, cos_dp)])

        dpfinal[index] = np.arctan2(
                sin_avg[index], cos_avg[index])*u.rad.to(u.deg)/file_factor

        if do_plot:
            plot(index, dp, dpfinal, wavelength,
                 "Differential phase", lim)
    return dpfinal


def calibrate_visibilities(
        va: List[np.ndarray], wavelength: np.ndarray,
        file_factor: int, do_plot: Optional[bool] = True) -> np.ndarray:
    """Calibrates the visibilities or correlated fluxes."""
    va = [compact_exposures(elem, compact_factor=6) for elem in va]
    vafinal = np.zeros((6, va[0].shape[1]))

    if do_plot:
        plt.figure(53)

    for index in range(6):
        if len(va) == 4:
            vafinal[index] = sum([elem[DP_INDEX[index, i]]
                                  for i, elem in enumerate(va)])
        else:
            indices = [0, 3]
            vafinal[index] = sum([elem[DP_INDEX[index, i]]
                                  for i, elem in zip(indices, va)])

        if do_plot:
            plot(index, va, vafinal, wavelength,
                 "VISAMP", [None, None])
    return vafinal/file_factor


def calibrate_squared_visibilities(
        v2: List[np.ndarray], wavelength: np.ndarray,
        file_factor: int, do_plot: Optional[bool] = True) -> np.ndarray:
    """Calibrates the squared visibilities."""
    v2 = [compact_exposures(elem, compact_factor=6) for elem in v2]
    v2final = np.zeros((6, v2[0].shape[1]))

    if do_plot:
        plt.figure(50)

    for index in range(6):
        if len(v2) == 4:
            v2final[index] = sum([elem[DP_INDEX[index, i]]
                                  for i, elem in enumerate(v2)])
        else:
            indices = [0, 3]
            v2final[index] = sum([elem[DP_INDEX[index, i]]
                                  for i, elem in zip(indices, v2)])

        if do_plot:
            plot(index, v2, v2final, wavelength,
                 "Squared visibility", [0, 1])
    return v2final/file_factor


def calib_BCD(iifile: Path, iofile: Path,
              oifile: Path, oofile: Path,
              outputfile: Optional[Path] = None,
              lim: Optional[int] = 180,
              do_plot: Optional[bool] = True) -> None:
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

    files = list(map(Path, [iifile, oofile]))
    if iofile not in [None, "", " "]:
        files.insert(2, iofile)
    if oifile not in [None, "", " "]:
        files.insert(2, oifile)

    t3p, dp, va, v2, wl, vis_type = zip(*map(read_fits, files))
    file_factor = sum([int(elem.shape[0]/6) for elem in va])

    closfinal = calibrate_closure_phases(t3p, wl[0], file_factor, lim, do_plot)
    dpfinal = calibrate_differential_phase(dp, wl[0], file_factor, lim, do_plot)
    vafinal = calibrate_visibilities(va, wl[0], file_factor, do_plot)
    v2final = calibrate_squared_visibilities(v2, wl[0], file_factor, do_plot)

    with fits.open(outputfile, mode="update") as hdul:
        hdul["OI_T3"].data = hdul["OI_T3"].data[0:4]
        hdul["OI_T3"].data["T3PHI"] = closfinal[0:4]
        hdul["OI_VIS"].data = hdul["OI_VIS"].data[0:6]
        hdul["OI_VIS"].data["VISPHI"] = dpfinal[0:6]
        hdul["OI_VIS"].data["VISAMP"] = vafinal
        hdul["OI_VIS2"].data = hdul["OI_VIS2"].data[0:6]
        hdul["OI_VIS2"].data["VIS2DATA"] = v2final

        if "correlated" == vis_type[0]:
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

    if do_plot:
        plt.show()
