import keyring
import time
from datetime import datetime, timedelta
from functools import wraps
from pathlib import Path
from typing import Union, Optional, Callable, Tuple, List

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import EarthLocation
from astropy.io import fits
from astropy.table import Table, vstack
from astroquery.eso import Eso
from scipy.interpolate import CubicSpline
from matplotlib import colormaps as mcm
from matplotlib.colors import ListedColormap

from .options import OPTIONS


__all__ = ["cprint", "capitalise_to_index", "print_execution_time",
           "get_execution_modes", "split_fits", "get_fits_by_tag",
           "check_if_target", "get_path_descriptor"]


class HeaderNotFoundWarning(Warning):
    """This indicates that a header of a (.fits)-file was not found."""
    pass


def convert_style_to_colormap(style: str) -> ListedColormap:
    """Converts a style into a colormap."""
    plt.style.use(style)
    colormap = ListedColormap(
            plt.rcParams["axes.prop_cycle"].by_key()["color"])
    plt.style.use("default")
    return colormap


def get_colormap(colormap: str) -> ListedColormap:
    """Gets the colormap as the matplotlib colormaps or styles."""
    try:
        return mcm.get_cmap(colormap)
    except ValueError:
        return convert_style_to_colormap(colormap)


def get_colorlist(colormap: str,
                  ncolors: Optional[int] = OPTIONS.color.number) -> List[str]:
    """Gets the colormap as a list from the matplotlib colormaps."""
    return [get_colormap(colormap)(i) for i in range(ncolors)]


# TODO: Get a better error representation for the flux.
# TODO: Implement smoothing for the flux to the instrument
def get_flux_data_from_flux_file(
        flux_file: Path, wavelength_axis: np.ndarray,
        error_percentage: float
        ) -> Tuple[u.Quantity[u.Jy], u.Quantity[u.Jy]]:
    """Reads the flux data from the flux file and then interpolates it
    to the wavelength solution used by MATISSE.

    Parameters
    ----------
    flux_file : pathlib.Path
        The flux file to be read.
    wavelength_axis : numpy.ndarray
        The wavelength axis to interpolate the flux file to.
    error_percentage : float
        The percentage of the flux error to be used in the interpolation

    Returns
    -------
    flux : astropy.units.Jy
        The, to MATISSE's wavelength solution interpolated, flux fetched
        from the flux file.
    flux_error : astropy.units.Jy
        The, to MATISSE's wavelength solution interpolated, flux error
        fetched from the flux file.
    """
    flux_data = Table.read(flux_file, names=["wl", "flux"], format="ascii")
    cubic_spline = CubicSpline(flux_data["wl"], flux_data["flux"])
    interpolated_flux = (cubic_spline(wavelength_axis)).ravel()
    return interpolated_flux, interpolated_flux*error_percentage


def unwrap_phases(phase: Union[float, np.ndarray],
                  error: Optional[Union[float, np.ndarray]] = None,
                  period: Optional[int] = 360
                  ) -> Union[Tuple[np.ndarray, np.ndarray], np.ndarray]:
    """Unwraps both the phases and the errors of the input arrays.

    Parameters
    ----------
    phase : float or numpy.ndarray
    error : float or numpy.ndarray, optional
    period : int, optional

    Returns
    -------
    phase : numpy.ndarray
    error : numpy.ndarray, optional
    """
    if error is not None:
        return map(lambda x: np.unwrap(x, period=period), [phase, error])
    return np.unwrap(phase, period=period)


def flip_phases(fits_file: Path) -> None:
    """Flips the phase of the t3phi and writes it in the header.

    Required for old versions (<=2.0.0) of the MATISSE data reduction pipeline
    which had flipped the N-band phase only by 180 degrees.
    """
    with fits.open(fits_file, "update") as hdul:
        header = hdul["oi_t3"].header
        if "PFLIP" in header:
            return
        hdul["oi_t3"].data["t3phi"] = -hdul["oi_t3"].data["t3phi"]
        header["PFLIP"] = True, "Phase flipped by 180 degrees"
        hdul.flush()


def cprint(message: str, color: Optional[str] = None) -> None:
    """Makes use of ascii-codes to print messages in color.

    Parameters
    ----------
    message : str
        The message to be printed.
    color : str, optional
        The name of the color to be used.
    """
    color_dict = {"r": ["\033[91m", "\033[00m"], "g": ["\033[92m", "\033[00m"],
                  "y": ["\033[93m", "\033[00m"], "lp": ["\033[94m", "\033[00m"],
                  "p": ["\033[95m", "\033[00m"], "cy": ["\033[96m", "\033[00m"],
                  "lg": ["\033[97m", "\033[00m"]}

    if color is not None:
        color_code = color_dict[color]
        print(color_code[0] + message + color_code[1])
    else:
        print(message)


def capitalise_to_index(string: str, index: int):
    """Capitalises a string until a certain index."""
    if index == len(string):
        return string.capitalize()
    return string[:index].capitalize() + string[index:].lower()


def print_execution_time(func: Callable):
    """Prints the execution time of the input function."""
    @wraps(func)
    def inner(*args, **kwargs):
        overall_start_time = time.perf_counter()
        result = func(*args, **kwargs)
        execution_time = time.perf_counter()-overall_start_time
        cprint(f"{'':-^50}", "lg")
        cprint(f"[INFO]: Executed in {timedelta(seconds=execution_time)}"
               " hh:mm:ss", "lg")
        cprint(f"{'':-^50}", "lg")
        return result
    return inner


def get_execution_modes(mode: Optional[str] = None,
                        band: Optional[str] = None) -> Tuple[List[str], List[str]]:
    """Determines the mode- and band configurations used by the users input. Returns
    either one or two lists depending on the input.

    Parameters
    ----------
    mode : str, optional
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'.
    band : str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'.

    Returns
    -------
    modes : list of str
        A list of the modes to be enhanced.
    bands : list of str
        A list of the bands to be enhanced.
    """
    modes, bands = [], []
    if mode is not None:
        if mode not in ["both", "coherent", "incoherent"]:
            raise IOError(f"No mode named '{mode}' exists!")
        modes = ["coherent", "incoherent"] if mode == "both" else [mode]

    if band is not None:
        if band not in ["both", "lband", "nband"]:
            raise IOError(f"No band named '{band}' exists!")
        bands = ["lband", "nband"] if band == "both" else [band]
    return modes, bands


def split_fits(directory: Path, tag: str) -> Tuple[List[Path], Optional[List[Path]]]:
    """Searches a folder for a tag and if files are found it returns the non-chopped
    and chopped (.fits)-files. If there are only non-chopped (.fits)-files it will return
    'None' for the chopped-files.

    Parameters
    ----------
    directory : pathlib.Path
        The directory to be searched in.
    tag : str
        The tag that is contained in the file names.

    Returns
    ----------
    unchopped_fits : list of pathlib.Path
        A list of Paths that are the chopped (.fits)-files.
    chopped_fits : list of pathlib.Path, optional
        A list of Paths that are the unchopped (.fits)-files.
    """
    fits_files = [fits_file for fits_file in get_fits_by_tag(directory, tag)]
    unchopped_fits, chopped_fits = [], []
    for fits_file in fits_files:
        index = int(fits_file.name[-6])
        if 0 <= index <= 4:
            unchopped_fits.append(fits_file)
        elif 5 <= index <= 6:
            chopped_fits.append(fits_file)
    return unchopped_fits, chopped_fits if chopped_fits else None


# TODO: Write function in such a way that it gets the name of the files from the headers
# and then checks if they are correct
def get_fits_by_tag(directory: Path, tag: str) -> List[Path]:
    """Searches a folder for a tag and returns the (.fits)-files matching it.

    Parameters
    ----------
    directory : pathlib.Path
        The directory to be searched in.
    tag : str
        The tag that is contained in the file names.

    Returns
    ----------
    files : list of pathlib.Path
        A list of files that contain the tag in their names
    """
    return sorted(directory.glob(f"*{tag}*.fits"), key=lambda x: x.name[-8:])


def check_if_target(target_dir: Path) -> bool:
    """Checks if the given directory contains 'TARGET_RAW_INT'-files.

    Parameters
    ----------
    target_dir : pathlib.Path
        The directory that is to be checked for 'TARGET_RAW_INT' files.

    Returns
    -------
    contains_target : bool
    """
    return True if target_dir.glob("TARGET_RAW_INT*") else False


def get_path_descriptor(root_dir: Path, descriptor: Path,
                        tar_dir: Path, cal_dir: Path) -> Path:
    """Assembles the names for the new directories that will contain the
    calibrated files and returns the 'output_dir'.

    Parameters
    ----------
    root_dir : pathlib.Path
        The root directory of the PRODUCT.
    descriptor : pathlib.Path
        The desired name 'TAR-CAL' directory.
    tar_dir : pathlib.Path
        The target directory.
    cal_dir : pathlib.Path
        The calibrator directory.

    Returns
    -------
    output_dir : pathlib.Path
        The path for the output directory.
    """
    mode_and_band = str(tar_dir.parents[1]).split("/")[-2:]
    dir_name, time_stamp_sci, detector = str(tar_dir.parent).split(".")[:-1]
    dir_name = dir_name.split('/')[-1].replace("raw", "cal")
    time_stamp_cal = str(cal_dir.parent).split('.')[-3]
    new_dir_name = '.'.join([descriptor, dir_name,
                             time_stamp_sci, detector, time_stamp_cal, "rb"])
    return root_dir / "calib" / mode_and_band[0] / new_dir_name


# TODO: Reimplement both of the following functions in a better way
def transform_uv_points(
        baselines: List[float], hour_angle: np.ndarray,
        latitude: u.rad, declination: u.rad) -> Tuple[np.ndarray, np.ndarray]:
    """Calculates the earth rotation (synthesis) for the uv-point
    corresponding to the baselines for the input hour angle(s)

    Parameters
    -----------
    baselines : list of float
        The baselines in the following order: Baselines east, -north, -longest.
    hour_angle : numpy.ndarray of float
    latitude : astropy.units.rad
        The latitude of the site
    declination : astropy.units.rad

    Returns
    -------
    u_coords : numpy.ndarray
    v_coords : numpy.ndarray
    """
    baseline_east, baseline_north, baseline_longest = baselines

    u_coords = baseline_east * np.cos(hour_angle) - baseline_north * np.sin(latitude)\
        * np.sin(hour_angle) + baseline_longest * np.cos(latitude) * np.sin(hour_angle)
    v_coords = baseline_east * np.sin(declination) * np.sin(hour_angle)\
        + baseline_north * (np.sin(latitude) * np.sin(declination) * np.cos(hour_angle)
                            + np.cos(latitude) * np.cos(declination)) - baseline_longest * \
        (np.cos(latitude) * np.sin(declination) * np.cos(hour_angle)
         - np.sin(latitude) * np.cos(declination))
    return u_coords, v_coords


def calculate_uv_tracks(baselines: List[np.ndarray],
                        declination: float, airmass_lim: float):
    """Calculates the tracks that uv-coordinates create from the earth rotation
    synthesis at the latitude of paranal.

    Parameters
    ----------
    baselines : list of numpy.ndarray
        The baselines in the following order: Baselines east, -north, -longest.
    declination : float
        The declination of the target.
    airmass_lim : float
        The airmass limit of the target.
    """
    latitude_paranal = EarthLocation.of_site(
        "paranal").geodetic.lat.to(u.rad)
    hamax = np.arccos(abs((1./airmass_lim-np.sin(latitude_paranal)
                           * np.sin(declination))/(np.cos(latitude_paranal)
                                                   * np.cos(declination))))
    return transform_uv_points(baselines, np.linspace(-hamax, hamax, 1000),
                               latitude_paranal, declination)


def query_archive(
        user_name: str, instrument: str,
        target: str, nights: List[str],
        columns: Optional[List[str]] = [""],
        remove_password: Optional[bool] = False,
        row_limit: Optional[int] = 1500,
        help: Optional[bool] = False,
        download: Optional[bool] = False,
        save_dir: Optional[str] = None
        ) -> None:
    """Queries the ESO's raw data archive.

    Parameters
    ----------
    user_name : str
        The user name for the ESO archive.
    instrument : str
        The instrument to query the archive for.
    target : str
        The target to query the archive for.
    nights : list of str
        The nights to be queried.
    row_limit : int, optional
        The maximum number of rows to return.
        If -1, all rows are returned (takes quite some time).
    remove_password : bool, optional
        If True, the password will be removed from the keyring.
    help : bool, optional
        If True, the help message will be printed.
    download : bool, optional
        If True, the data will be downloaded.
    save_dir : str, optional
        The directory to save the data to. The data will be sorted in
        the folders pertaining to their array configuration and night.
    """
    eso = Eso()
    eso.ROW_LIMIT = row_limit
    eso.login(user_name, store_password=True)

    if remove_password:
        keyring.delete_password("astroquery:www.eso.org", user_name)

    for night in nights:
        next_night = (datetime.strptime(night, "%Y-%m-%d")+timedelta(days=1)).strftime("%Y-%m-%d")
        column_filters = {"target": target, "stime": night, "etime": next_night}
        sci = eso.query_instrument(instrument, column_filters=column_filters,
                                   columns=columns, help=help)
        sci_headers = eso.get_headers(sci['DP.ID'])
        sci_starts = np.unique(sci_headers["HIERARCH ESO TPL START"]).data
        sci_name = sci_headers["HIERARCH ESO OBS TARG NAME"][0]
        array_configuration = sci_headers["HIERARCH ESO ISS BASELINE"][0]

        del column_filters["target"]
        night = eso.query_instrument(instrument, column_filters=column_filters,
                                     columns=columns, help=help)
        night_headers = eso.get_headers(night['DP.ID'])
        tpl_starts = np.unique(night_headers["HIERARCH ESO TPL START"]).data
        ind = int(np.where(tpl_starts == sci_starts)[0].squeeze())

        calibrators = []
        for i in [ind-1, ind+1]:
            obs = night_headers[night_headers["HIERARCH ESO TPL START"] == tpl_starts[i]]
            dpr_type = obs["HIERARCH ESO DPR CATG"][0]
            if dpr_type == "SCIENCE":
                if obs["HIERARCH ESO OBS TARG NAME"][0] == sci_name:
                    print("[WARNING]: Sequential observation of the same SCI.")
                continue
            else:
                ob_name = obs["HIERARCH ESO OBS NAME"][0]
                calibrators.append(obs)
                if sci_name.replace("_", "").lower() in ob_name.replace("_", "").lower():
                    if "LN" in ob_name.upper():
                        break
    
        if download:
            download_files = vstack([sci_headers, *calibrators])["DP.ID"]
            destination = save_dir / array_configuration / night

            eso.retrieve_data(
                    download_files, destination=destination,
                    with_calib="raw", continuation=True,
                    request_all_objects=True, unzip=True)

