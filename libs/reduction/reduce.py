import os
import time
import datetime
import shutil
from pathlib import Path
from typing import Set, Optional

import numpy as np
import astropy.units as u
from astropy.table import Table
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from mat_tools import mat_autoPipeline as mp

from ..utils.plot import Plotter
from ..utils.readout import ReadoutFits
from ..utils.tools import cprint, get_fits_by_tag


DATA_DIR = Path(__file__).parent.parent.parent / "data"
CATALOG_DIR = DATA_DIR / "catalogues"

JSDC_V2_CATALOG = Vizier(catalog="II/346/jsdc_v2")
JSDC_CATALOG = CATALOG_DIR / "jsdc_v2_catalog_20170303.fits"
ADDITIONAL_CATALOG = CATALOG_DIR / "supplementary_catalog_202207.fits"

SPECTRAL_BINNING = {"low": [5, 7], "high_ut": [5, 38], "high_at": [5, 98]}


# TODO: Maybe there is a smarter way than globbing twice? Not that important however...
def get_tpl_match(raw_dir: Path, tpl_start: str) -> Path:
    """Gets a singular (.fits)-file matching the 'tpl_start', i.e., the starting time of
    the observation

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
        if tpl_start == ReadoutFits(fits_file).tpl_start:
            return fits_file
    raise FileNotFoundError("No file with matching tpl_start exists!")


def get_tpl_starts(raw_dir: Path) -> Set[str]:
    """Iterates through all files and gets their 'tpl_start', i.e, the starting time of
    the individual observations

    Parameters
    ----------
    raw_dir: Path
        The directory containing the raw-files
    """
    return set([ReadoutFits(fits_file).tpl_start for fits_file in raw_dir.glob("*.fits")])


def in_catalog(readout: ReadoutFits, radius: u.arcsec, catalog: Path) -> Optional[Path]:
    """Checks if calibrator is in the given catalog.

    Parameters
    ----------
    readout: ReadoutFits
        A class wrapping a (.fits)-file, that reads out its information quickly
    radius: u.arcsec
        The radius in which targets are queried from the catalog
    catlog: Path
        The catalog which is to be queried from the catalog

    Returns
    -------
    catalog: Path | None
        The catalog in which the object has been found in
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


def get_catalog_match(readout: ReadoutFits, radius: u.arcsec = 20*u.arcsec):
    """Checks if the calibrator is in the 'jsdc_v2'-catalog and if not then searches the
    local, supplementary calibrator databases

    Parameters
    ----------
    readout: ReadoutFits
        A class wrapping a (.fits)-file, that reads out its information quickly
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
            cprint(f"Calibrator '{match[0]['Name'][0]}' found in JSDC v2 catalog!",
                   "y")
        return JSDC_CATALOG
    return in_catalog(readout.coords, radius=radius, catalog=ADDITIONAL_CATALOG)


# TODO: Determine the spectral binning by automatic readout in the future
def set_script_arguments(corr_flux: bool, array: str,
                         resolution: Optional[str] = "low") -> Tuple[str]:
    """Sets the arguments that are then passed to the 'mat_autoPipeline.py'
    script

    Parameters
    ----------
    corr_flux: bool
        This specifies if the flux is to be reduced or not
    array: str
        The array configuration that was used for the observation
    resolution: str, optional
        This determines the spectral binning. Input can be "low" for
        low-resolution in both bands, "high_ut" for low-resolution in L-band
        and high-resolution in N-band for the UTs and the same for "high_at" for
        the ATs

    Returns
    -------
    lband_params: str
        The arguments passed to the MATISSE-pipline for the L-band
    nband_params: str
        The arguments passed to the MATISSE-pipline for the N-band
    """
    if resolution == "high":
        resolution = f"{resolution}_{array.lower()[:-1]}"
    bin_L, bin_N = SPECTRAL_BINNING[resolution]

    # NOTE: Jozsef uses TEL 3 here, but Jacob 2? What is the meaning of this
    # -> Read up on it. Already asked! Awaiting response
    compensate = '/compensate="[pb,rb,nl,if,bp,od]"'
    tel = "/replaceTel=3" if array == "ATs" else "/replaceTel=0"
    coh_L = "/corrFlux=TRUE/useOpdMod=FALSE/coherentAlgo=2" if corr_flux else ""
    coh_N = "/corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2" if corr_flux else ""
    return f"{coh_L}{compensate}/spectralBinning={bin_L}",\
            f"{tel}{coh_N}/spectralBinning={bin_N}"


# TODO: Read out the array configurations from the (.fits)-files
def reduce_mode_and_band(raw_dir: Path, calib_dir: Path, res_dir: Path,
                         array: str, mode: bool, band: str, tpl_start: str,
                         resolution: Optional[str] = "low") -> None:
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
        The directory containing the calibration files
    res_dir: Path
        The directory to contain to reduced data
    array: str
        The array configuration that was used for the observation. Either "AT" or "UT"
    mode: bool
        The mode in which the reduction is to be done, either "incoherent" if "False" or
        "coherent" if "True"
    band: str
        The band for which the reduction is to be done, either "lband" or "nband"
    tpl_star: str
        The starting time of observations
    """
    start_time = time.perf_counter()
    mode_and_band_dir = res_dir / mode / band
    param_L, param_N = set_script_arguments(mode, array, resolution)
    skip_L = True if band == "nband" else False
    skip_N = not skip_L

    if not mode_and_band_dir.exists():
        mode_and_band_dir.mkdir(parents=True)

    readout = ReadoutFits(get_tpl_match(raw_dir, tpl_start))
    if readout.is_calibrator():
        cprint(f"Calibrator '{readout.name}' detected!"
               f" Checking for catalog...", "g")
        catalog_path = get_catalog_match(readout)
        if catalog_path is not None:
            cprint(f"Moving catalog to {calib_dir.name}...", "g")
            shutil.copy(catalog_path, calib_dir / catalog_path.name)
    else:
        for catalog in calib_dir.glob("*catalog*"):
            os.remove(catalog)
        cprint(f"Science target '{readout.name}' detected! Removing catalogs...", "y")

    mp.mat_autoPipeline(dirRaw=str(raw_dir), dirResult=str(res_dir),
                        dirCalib=str(calib_dir), tplstartsel=tpl_start,
                        nbCore=6, resol='', paramL=param_L, paramN=param_N,
                        overwrite=0, maxIter=1, skipL=skip_L, skipN=skip_N)

    try:
        rb_folders = res_dir.glob("Iter1/*.rb")
        for folder in rb_folders:
            cprint(f"Moving folder {folder.name}...", "g")
            if (mode_and_band_dir / folder.name).exists():
                shutil.rmtree(mode_and_band_dir / folder.name)
            shutil.move(folder, mode_and_band_dir)

        if rb_folders:
            cprint(f"Moving folders to directory {mode_and_band_dir.name}...",
                   "g")

    # TODO: Make logger here
    except Exception:
        cprint(f"Moving of files to {mode_and_band_dir.name} failed!", "r")

    for folder in mode_and_band_dir.glob("*.rb"):
        cprint(f"Plotting files of folder {folder.name}...", "g")
        for fits_file in get_fits_by_tag(folder, "RAW_INT"):
            plot_fits = Plotter([fits_file],
                                save_path=(mode_and_band_dir / folder.name))
            plot_fits.add_cphase().add_vis().plot(save=True)

    cprint(f"{'':-^50}", "lg")
    cprint(f"Executed the {mode} reduction for the {band} in"
           f" {datetime.timedelta(seconds=(time.perf_counter()-start_time))}"
           "hh:mm:ss",
           "lg")
    cprint(f"{'':-^50}", "lg")


# TODO: Change the way the folders are structured and input for ease-of use
def reduce(root_dir: Path, instrument_dir: Path,
           target_dir: Path, array: str,
           resolution: Optional[str] = "low"):
    """Runs the pipeline for the data reduction

    Parameters
    ----------
    raw_dir: Path
        The directory containing the raw observation files
    array: str
        The array configuration that was used for the observation
    """
    overall_start_time = time.perf_counter()
    raw_dir = Path(root_dir, instrument_dir, "raw", target_dir).resolve()

    calib_dir = raw_dir / "calib_files"
    # TODO: Catch edge case where this does not move the files if calib_dir
    # exists already!
    if not calib_dir.exists():
        calib_dir.mkdir(parents=True)
        cprint("Moving calibration files into 'calib_files' folders...", "g")
        calibration_files = raw_dir.glob("M.*")
        for calibration_file in calibration_files:
            shutil.move(calibration_file, calib_dir / calibration_file.name)

    res_dir = Path(root_dir, instrument_dir, "products", target_dir).resolve()

    if not res_dir.exists():
        res_dir.mkdir(parents=True)

    # TODO: Add in the option to not remove old reduction and make new one take an
    # addional tag after its name
    try:
        folders = res_dir.glob("*")
        for folder in folders:
            if folder.exists():
                shutil.rmtree(folder)
        cprint("Cleaning up old reduction...", "y")

    # TODO: Make logger here
    except Exception:
        cprint("Cleaning up failed!", "y")

    for tpl_start in sorted(list(get_tpl_starts(raw_dir))):
        cprint(f"{'':-^50}", "lg")
        cprint(f"Reducing data of tpl_start: {tpl_start}", "g")
        cprint(f"{'':-^50}", "lg")
        for mode in ["coherent", "incoherent"]:
            cprint(f"Processing {mode} reduction...", "lp")
            cprint(f"{'':-^50}", "lg")
            for band in ["lband", "nband"]:
                reduce_mode_and_band(raw_dir, calib_dir, res_dir, array,
                                     resolution=resolution, mode=mode, band=band,
                                     tpl_start=tpl_start)

    execution_time = time.perf_counter()-overall_start_time
    cprint(f"{datetime.timedelta(seconds=execution_time)} hh:mm:ss", "lg")


if __name__ == "__main__":
    root_dir = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/")
    target_dir, date_dir = "matisse/GTO/hd142666", "UTs/20190514"
    reduce(root_dir, target_dir, date_dir, "UTs")
