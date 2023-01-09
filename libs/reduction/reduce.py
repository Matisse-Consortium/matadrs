import os
import time
import datetime
import shutil
from pathlib import Path
from typing import Optional

import numpy as np
import astropy.units as u
from astropy.table import Table
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from mat_tools import mat_autoPipeline as mp

# TODO: Find way to make this into a complete module -> More pythonic!
from utils import cprint
from readout import ReadoutFits


DATA_DIR = Path(__file__).parent.parent.parent / "data"
CATALOG_DIR = DATA_DIR / "catalogues"

JSDC_V2_CATALOG = Vizier(catalog="II/346/jsdc_v2")
JSDC_CATALOG = CATALOG_DIR / "jsdc_v2_catalog_20170303.fits"
ADDITIONAL_CATALOG = CATALOG_DIR / "supplementary_catalog_202207.fits"

SPECTRAL_BINNING = {"low": [5, 7], "high_ut": [5, 38], "high_at": [5, 98]}


# TODO: Maybe there is a smarter way than globbing twice? Not that important however...
def get_tpl_match(raw_dir: Path, tpl_start: str):
    """Gets a (.fits)-file matching the tpl start

    Parameters
    ----------
    raw_dir: Path
    tpl_start: str

    Returns
    -------
    fits_file: Path
    """
    for fits_file in raw_dir.glob("*.fits"):
        if tpl_start == ReadoutFits(fits_file).tpl_start:
            return fits_file
    raise FileNotFoundError("No file with matching tpl_start exists!")


def get_tpl_starts(raw_dir: Path):
    """Iterates through all files and gets their tpl start times"""
    return set([ReadoutFits(fits_file).tpl_start for fits_file in raw_dir.glob("*.fits")])


def in_catalog(coords_calibrator: SkyCoord, radius: u.arcsec, catalog: Path):
    """Checks if calibrator is in catalog. Returns catalog if True, else None

    Parameters
    ----------
    coords_calibrator: SkyCoord
    radius: u.arcsec
    catalog_fits: Path

    Returns
    -------
    catalog: Path | None
    """
    table = Table().read(catalog)
    coords_catalog = SkyCoord(table["RAJ2000"], table["DEJ2000"],
                              unit=(u.hourangle, u.deg), frame="icrs")
    separation = coords_calibrator.separation(coords_catalog)
    if separation[np.nanargmin(separation)] < radius.to(u.deg):
        cprint("Calibrator found in supplementary catalog!", "g")
        return catalog
    cprint("Calibrator not found in any catalogs!"
           " No TF2 and 'mat_cal_oifits' will fail", "r")
    return None


def get_catalog_match(fits_file: Path, match_radius: u.arcsec = 20*u.arcsec):
    """Checks if the calibrator is in the 'jsdc_v2'-catalog and if not then searches the
    local calibrator databases

    Parameters
    ----------
    fits_file: Path
        The file that is to be reduced and checked if calibrator or not
    match_radius: u.arcsec
        The radius in which to search the catalouge

    Returns
    -------
    catalog: Path | None
    """
    coords_calibrator = ReadoutFits(fits_file).coords
    if JSDC_V2_CATALOG.query_region(coords_calibrator, radius=match_radius):
        cprint("Calibrator found in JSDC v2 catalog!", "g")
        return JSDC_CATALOG
    return in_catalog(coords_calibrator, radius=match_radius, catalog=ADDITIONAL_CATALOG)


def set_script_arguments(corr_flux: bool, array: str,
                         resolution: Optional[str] = "low") -> str:
    """Sets the arguments that are then passed to the 'mat_autoPipeline.py'
    script

    Parameters
    ----------
    corr_flux: bool
        This specifies if the flux is to be reduced or not
    array: str
        The array configuration that was used for the observation
    resolution: str
        This determines the spectral binning. Input can be "low" for
        low-resolution in both bands, "high_ut" for low-resolution in L-band
        and high-resolution in N-band for the UTs and the same for "high_at" for
        the ATs

    Returns
    -------
    str
        A string that contains the arguments, which are passed to the
        MATISSE-pipline
    """
    bin_L, bin_N = SPECTRAL_BINNING[resolution]
    # NOTE: Jozsef uses TEL 3 here, but Jacob 2? What is the meaning of this
    # -> Read up on it. Already asked! Awaiting response
    compensate = "/compensate=[pb,rb,nl,if,bp,od]"
    tel = "/replaceTel=3" if array == "ATs" else "/replaceTel=0"
    coh_L  = f"/corrFlux=TRUE/useOpdMod=FALSE/coherentAlgo=2"\
            if corr_flux else ""
    coh_N = f"/corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2" if corr_flux else ""
    paramL_lst = f"{coh_L}{compensate}/spectralBinning={bin_L}"
    paramN_lst = f"{tel}{coh_N}/spectralBinning={bin_N}"
    return paramL_lst, paramN_lst


def reduce_mode_and_band(raw_dir: Path, calib_dir: Path, res_dir: Path,
                         array: str, mode: bool, band: str, tpl_start: str,
                         resolution: Optional[str] = "low") -> None:
    """Reduces either the lband or the nband data for either the "coherent" or
    "incoherent" setting for a single iteration/epoch.

    Removes the old (.sof)-files to ensure proper reduction and then creates the needed
    folders in the "res_dir"-directory. After this, it starts the reduction with the
    specified settings

    Parameters
    ----------
    raw_dir: Path
        The path containing the raw observation files
    calib_dir: Path
        The path containing the calibration files
    res_dir: Path
        The path to contain to reduced data
    array: str
        The array configuration that was used for the observation. Either "AT" or "UT"
    mode: bool
        The mode in which the reduction is to be done, either "incoherent" if "False" or
        "coherent" if "True"
    band: str
        The band for which the reduction is to be done, either "lband" or "nband"
    tpl_star: str
        The starting time of target's observations
    """
    start_time = time.perf_counter()
    mode_and_band_dir = res_dir / mode / band
    param_L, param_N = set_script_arguments(mode, array, resolution)
    skip_L = True if band == "nband" else False
    skip_N = not skip_L

    if not mode_and_band_dir.exists():
        mode_and_band_dir.mkdir(parents=True)

    # TODO: Readout is called twice here, make this easier!!!
    fits_file = get_tpl_match(raw_dir, tpl_start)
    if ReadoutFits(fits_file).is_calibrator():
        cprint("Calibrator detected!"
               f" Moving corresponding catalog to {calib_dir.name}...", "y")
        catalog_path = get_catalog_match(fits_file)
        if catalog_path is not None:
            shutil.copy(catalog_path, calib_dir / catalog_path.name)
    else:
        for catalog in calib_dir.glob("*catalog*"):
            os.remove(catalog)
        cprint("Science target detected! Removing catalogs...", "y")

    mp.mat_autoPipeline(dirRaw=str(raw_dir), dirResult=str(res_dir),
                        dirCalib=str(calib_dir), tplstartsel=tpl_start,
                        nbCore=6, resol='', paramL=param_L, paramN=param_N,
                        overwrite=0, maxIter=1, skipL=skip_L, skipN=skip_N)

    try:
        rb_folders = res_dir.glob("Iter1/*.rb")
        for folder in rb_folders:
            if (mode_and_band_dir / folder.name).exists():
                shutil.rmtree(mode_and_band_dir / folder.name)
            shutil.move(folder, mode_and_band_dir)

        if rb_folders:
            cprint(f"Moving folders to directory {mode_and_band_dir.name}", "g")

    # TODO: Make logger here
    except Exception:
        cprint(f"Moving of files to {mode_and_band_dir.name} failed!", "y")

    cprint(f"{'':-^50}", "lg")
    cprint(f"Executed the {mode} reduction for the {band} in"
           f" {datetime.timedelta(seconds=(time.perf_counter()-start_time))}"
           "hh:mm:ss",
           "lg")
    cprint(f"{'':-^50}", "lg")


def reduce(root_dir: Path, instrument_dir: Path,
           target_dir: Path, array: str,
           resolution: Optional[str] = "low"):
    """Runs the pipeline for the data reduction

    Parameters
    ----------
    raw_dir: Path
        The path containing the raw observation files
    array: str
        The array configuration that was used for the observation

    See Also
    --------
    single_reduction()
    """
    overall_start_time = time.perf_counter()
    raw_dir = Path(root_dir, instrument_dir, "raw", target_dir).resolve()

    calib_dir = raw_dir / "calib_files"
    # TODO: Catch edge case where this does not move the files if calib_dir
    # exists already!
    if not calib_dir.exists():
        calib_dir.mkdir(parents=True)
        cprint("Moving calibration files into 'calib_files' folders!", "g")
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
            cprint(f"Processing {mode} reduction", "lp")
            cprint(f"{'':-^50}", "lg")
            for band in ["lband", "nband"]:
                reduce_mode_and_band(raw_dir, calib_dir, res_dir, array,
                                     resolution=resolution, mode=mode, band=band,
                                     tpl_start=tpl_start)

    execution_time = time.perf_counter()-overall_start_time
    cprint(f"{datetime.timedelta(seconds=execution_time)} hh:mm:ss", "lg")


if __name__ == "__main__":
    root_dir = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/")
    target_dir, date_dir = "matisse/GTO/hd142666", "UTs/20220420"
    reduce(root_dir, target_dir, date_dir, "UTs")
