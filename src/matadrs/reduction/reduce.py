import os
import shutil
from pathlib import Path
from typing import Set, Tuple, Union, Optional

import numpy as np
import astropy.units as u
from astropy.table import Table
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord

from ..mat_tools.mat_autoPipeline import mat_autoPipeline as mp
from ..mat_tools.libAutoPipeline import matisseType
from ..utils.plot import Plotter
from ..utils.readout import ReadoutFits
from ..utils.tools import cprint, print_execution_time,\
        get_execution_modes, get_fits_by_tag

__all__ = ["get_readout_for_tpl_match", "get_tpl_starts", "in_catalog",
           "get_catalog_match", "prepare_catalogs", "set_script_arguments",
           "cleanup_reduction", "reduce_mode_and_band", "prepare_reduction", "reduce"]

# CATALOG_DIR = pkgutil.get_data(__name__, "data/catalogues")
JSDC_V2_CATALOG = Vizier(catalog="II/346/jsdc_v2")
JSDC_CATALOG = CATALOG_DIR / "jsdc_v2_catalog_20170303.fits"
ADDITIONAL_CATALOG = CATALOG_DIR / "supplementary_catalog_202207.fits"

SPECTRAL_BINNING = {"low": [5, 7], "high_uts": [5, 38], "high_ats": [5, 98]}


def get_readout_for_tpl_match(raw_dir: Path, tpl_start: str) -> Path:
    """Gets the readout of a singular (.fits)-file matching the 'tpl_start', i.e.,
    the starting time of the observation

    Parameters
    ----------
    raw_dir: Path
        The directory containing the raw-files
    tpl_start: str
        The starting time of the observation

    Returns
    -------
    fits_file: Path
        A (.fits)-file matching the input
    """
    for fits_file in raw_dir.glob("*.fits"):
        readout = ReadoutFits(fits_file)
        if tpl_start == readout.tpl_start:
            return readout
    raise FileNotFoundError(f"No file with matching tpl_start: '{tpl_start}' exists!")


def get_tpl_starts(raw_dir: Path) -> Set[str]:
    """Iterates through all files and gets their 'tpl_start', i.e, the starting time of
    the individual observations

    Parameters
    ----------
    raw_dir: Path
        The directory containing the raw-files

    Returns
    -------
    tpl_starts: Set[str]
        The starting times of all the observations given by the raw-files
    """
    return set([ReadoutFits(fits_file).tpl_start for fits_file in raw_dir.glob("*.fits")])


def remove_old_catalogs(catalog: Path, calib_dir: Path) -> bool:
    """Checks if the latest catalog is already existing in the calibration directory and
    removes outdated iterations

    Parameters
    ----------
    calib_dir: Path
        The directory containing to the observation associated calibration files

    Returns
    -------
    outdated_catalogs_found: bool
        Returns 'True' if outdated catalog has been found
    """
    catalog_header = ReadoutFits(catalog).primary_header
    for fits_file in calib_dir.glob("*.fits"):
        header = ReadoutFits(fits_file).primary_header
        if (matisseType(header) == "JSDC_CAT") and ("DATE" in (header and catalog_header)):
            if header["DATE"] == catalog_header["DATE"]:
                outdated_catalogs_found = False
        else:
            outdated_catalogs_found = True
            os.remove(str(fits_file))
    return outdated_catalogs_found


def in_catalog(readout: ReadoutFits,
               radius: u.arcsec, catalog: Path) -> Optional[Path]:
    """Checks if calibrator is in the given supplementary catalog

    Parameters
    ----------
    readout: ReadoutFits
        A class wrapping a (.fits)-file, that reads out its information
    radius: u.arcsec
        The radius in which targets are queried from the catalog
    catlog: Path
        The catalog which is to be queried from the catalog

    Returns
    -------
    catalog: Path | None
        The catalog in which the object has been found in.
        Returns None if object not found
    """
    table = Table().read(catalog)
    coords_catalog = SkyCoord(table["RAJ2000"], table["DEJ2000"],
                              unit=(u.hourangle, u.deg), frame="icrs")
    separation = readout.coords.separation(coords_catalog)
    if separation[np.nanargmin(separation)] < radius.to(u.deg):
        cprint(f"Calibrator '{readout.name}' found in supplementary catalog!", "g")
        return catalog
    cprint(f"Calibrator '{readout.name}' not found in any catalogs!"
           " No TF2 and the 'mat_cal_oifits'-recipe will fail", "r")
    return None


def get_catalog_match(readout: ReadoutFits,
                      radius: u.arcsec = 20*u.arcsec) -> Union[Path, None]:
    """Checks if the given calibrator is contained in the 'jsdc_v2'-catalog. If otherwise
    searches the local, supplementary calibrator databases instead. If nothing is found
    returns None

    Parameters
    ----------
    readout: ReadoutFits
        A class wrapping a (.fits)-file, that reads out its information
    radius: u.arcsec
        The radius in which targets are queried from the catalog

    Returns
    -------
    catalog: Path | None
        The catalog in which the object has been found in
    """
    match = JSDC_V2_CATALOG.query_region(readout.coords, radius=radius)
    if match:
        if len(match[0]) > 0:
            cprint(f"Calibrator '{match[0]['Name'][0]}' found in JSDC v2 catalog!", "y")
        return JSDC_CATALOG
    return in_catalog(readout.coords, radius=radius, catalog=ADDITIONAL_CATALOG)


# FIXME: Somehow the Science target is also detected to be a calibrator, or not?
def prepare_catalogs(raw_dir: Path, calib_dir: Path, tpl_start: str) -> None:
    """Checks if the starting time given by 'tpl_start' corresponds to the observation of
    a science target/calibrator and removes/prepares the un/needed catalogs

    Parameters
    ----------
    raw_dir: Path
        The direcotry containing the raw observation files
    calib_dir: Path
        The directory containing to the observation associated calibration files
    tpl_start: str
        The starting time of the observation
    """
    readout = get_readout_for_tpl_match(raw_dir, tpl_start)
    if readout.is_calibrator():
        cprint(f"Calibrator '{readout.name}' detected!"
               f" Checking for catalog...", "g")
        catalog = get_catalog_match(readout)
        if catalog is not None:
            outdated_catalogs_found = remove_old_catalogs(calib_dir)
            if outdated_catalogs_found:
                cprint(f"Moving catalog to {calib_dir.name}...", "g")
                shutil.copy(catalog, calib_dir / catalog.name)
            else:
                cprint("Latest catalog already present!", "g")
    else:
        for catalog in calib_dir.glob("*catalog*"):
            os.remove(catalog)
        cprint(f"Science target '{readout.name}' detected! Removing catalogs...", "y")


def set_script_arguments(raw_dir: Path, mode: str, tpl_start: str) -> Tuple[str]:
    """Sets the arguments that are then passed to the 'mat_autoPipeline.py' script

    Parameters
    ----------
    raw_dir: Path
        The directory containing the raw observation files
    mode: str
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'
    tpl_start: str
        The starting time of the observation

    Returns
    -------
    lband_params: str
        The arguments passed to the MATISSE-pipline for the L-band
    nband_params: str
        The arguments passed to the MATISSE-pipline for the N-band
    """
    readout = get_readout_for_tpl_match(raw_dir, tpl_start)
    if readout.resolution == "high":
        resolution = f"{readout.resolution}_{readout.array_configuration}"
    tel = "/replaceTel=3" if readout.array_configuration == "ats" else "/replaceTel=0"
    coh_lband = coh_nband = ""
    if mode == "coherent":
        coh_lband = "/corrFlux=TRUE/useOpdMod=FALSE/coherentAlgo=2"
        coh_nband = "/corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2"
    compensate = '/compensate="[pb,rb,nl,if,bp,od]"'
    bin_lband, bin_nband = SPECTRAL_BINNING[resolution]
    return f"{coh_lband}{compensate}/spectralBinning={bin_lband}",\
            f"{tel}{coh_nband}/spectralBinning={bin_nband}"


def cleanup_reduction(product_dir: Path, mode: str,
                      band: str, overwrite: bool) -> None:
    """Moves the folders to their corresponding folders of structure '/mode/band' after
    the reduction has been finished and plots the (.fits)-files contained in them

    Parameters
    ----------
    mode: str, optional
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'
    band: str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'
    overwrite: bool, optional
        If 'True' overwrites present files from previous reduction
    """
    mode_and_band_dir = product_dir / mode / band
    if not mode_and_band_dir.exists():
        mode_and_band_dir.mkdir(parents=True)

    for reduced_folder in product_dir.glob("Iter1/*.rb"):
        cprint(f"Moving folder '{reduced_folder.name}'...", "g")
        if ((mode_and_band_dir / reduced_folder.name).exists()) and (not overwrite):
            cprint(f"Could not move {reduced_folder},"
                   " as directories from previous reduction exists!")
        else:
            shutil.move(reduced_folder, mode_and_band_dir)

    # TODO: Remove this for loop? Maybe after properly implementing plotting?
    for reduced_folder in mode_and_band_dir.glob("*.rb"):
        cprint(f"Plotting files of folder {reduced_folder.name}...", "g")
        for fits_file in get_fits_by_tag(reduced_folder, "RAW_INT"):
            plot_fits = Plotter([fits_file],
                                save_path=(mode_and_band_dir / reduced_folder.name))
            plot_fits.add_cphase().add_vis().plot(save=True)


@print_execution_time
def reduce_mode_and_band(raw_dir: Path, calib_dir: Path,
                         product_dir: Path, mode: bool,
                         band: str, tpl_start: str, overwrite: bool) -> None:
    """Reduces either the L- and/or the N-band data for either the 'coherent' and/or
    'incoherent' setting for a single iteration/epoch.

    Notes
    -----
    Removes the old (.sof)-files to ensure proper reduction and then creates folders in
    the 'res_dir'-directory. After this, it starts the reduction with the specified
    settings

    Parameters
    ----------
    raw_dir: Path
        The direcotry containing the raw observation files
    calib_dir: Path
        The directory containing to the observation associated calibration files
    product_dir: Path
        The directory to contain the reduced files
    mode: str
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'
    band: str
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'
    tpl_star: str
        The starting time of observations
    overwrite: bool, optional
        If 'True' overwrites present files from previous reduction
    """
    skip_L = True if band == "nband" else False
    skip_N = not skip_L

    # FIXME: Get array and resolution from the (.fits)-file
    param_L, param_N = set_script_arguments(mode, mode, tpl_start)
    mp.mat_autoPipeline(dirRaw=str(raw_dir), dirResult=str(product_dir),
                        dirCalib=str(calib_dir), tplstartsel=tpl_start,
                        nbCore=6, resol='', paramL=param_L, paramN=param_N,
                        overwrite=0, maxIter=1, skipL=skip_L, skipN=skip_N)

    cleanup_reduction(product_dir, mode, band, overwrite)
    cprint(f"{'':-^50}", "lg")


def prepare_reduction(raw_dir: Path, calib_dir: Path,
                      product_dir: Path, remove_previous: bool) -> None:
    """Prepares the reduction by removing removing old product files and sorting the raw
    files by associated calibrations and observations

    Parameters
    ----------
    raw_dir: Path
        The direcotry containing the raw observation files
    calib_dir: Path
        The directory containing to the observation associated calibration files
    product_dir: Path
        The directory to contain the reduced files
    """
    if not calib_dir.exists():
        calib_dir.mkdir(parents=True)

    cprint("Moving calibration files into 'calib_files' folders...", "g")
    calibration_files = raw_dir.glob("M.*")
    for calibration_file in calibration_files:
        shutil.move(calibration_file, calib_dir / calibration_file.name)

    if not product_dir.exists():
        product_dir.mkdir(parents=True)


@print_execution_time
def reduce(raw_dir: Path, product_dir: Path, mode: Optional[str] = "both",
           band: Optional[str] = "both", overwrite: Optional[bool] = False) -> None:
    """Runs the pipeline for the data reduction

    Parameters
    ----------
    raw_dir: Path
        The directory containing the raw observation files
    product_dir: Path
        The directory to contain the reduced files
    mode: str, optional
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'
    band: str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'
    overwrite: bool, optional
        If 'True' overwrites present files from previous reduction

    Notes
    -----
    The reduction is executed on 6 cores in multiprocessing
    """
    raw_dir = Path(raw_dir).resolve()
    calib_dir = raw_dir / "calib_files"
    product_dir = Path(product_dir / "reduced").resolve()
    # TODO: Remove the old files during the loop? -> To make it more versatile?
    prepare_reduction(raw_dir, calib_dir, product_dir)
    modes, bands = get_execution_modes(mode, band)

    for tpl_start in sorted(list(get_tpl_starts(raw_dir))):
        cprint(f"{'':-^50}", "lg")
        cprint(f"Reducing data of tpl_start: {tpl_start}", "g")
        cprint(f"{'':-^50}", "lg")
        for mode in modes:
            cprint(f"Processing {mode} reduction...", "lp")
            cprint(f"{'':-^50}", "lg")
            for band in bands:
                reduce_mode_and_band(raw_dir, calib_dir, product_dir,
                                     mode, band, tpl_start, overwrite)


if __name__ == "__main__":
    breakpoint()
