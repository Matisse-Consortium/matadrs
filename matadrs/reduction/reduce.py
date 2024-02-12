import os
import shutil
import warnings
from pathlib import Path
from types import SimpleNamespace
from typing import Optional, Tuple, List, Set

import astropy.units as u
import numpy as np
import pkg_resources
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
from astroquery.vizier import Vizier
from tqdm import tqdm

from ..mat_tools.libAutoPipeline import matisseType
from ..mat_tools.mat_autoPipeline import mat_autoPipeline
from ..utils.plot import Plotter, plot_data_quality
from ..utils.readout import ReadoutFits
from ..utils.tools import HeaderNotFoundWarning, cprint, \
        print_execution_time, get_execution_modes, get_fits_by_tag


# NOTE: Remove the headers warnings as raw files are non-oifits
warnings.filterwarnings("ignore", category=HeaderNotFoundWarning)

__all__ = ["set_script_arguments", "cleanup_reduction",
           "reduce_mode_and_band", "prepare_reduction"]


CATALOG_DIR = Path(
    pkg_resources.resource_filename("matadrs", "data/catalogues"))
JSDC_V2_CATALOG = Vizier(catalog="II/346/jsdc_v2")
JSDC_CATALOG = CATALOG_DIR / "jsdc_v2_catalog_20170303.fits"
ADDITIONAL_CATALOG = CATALOG_DIR / "supplementary_catalog_202207.fits"

SPECTRAL_BINNING = {"low": [5, 7], "high_uts": [5, 38], "high_ats": [5, 98]}
CALIBRATION_IDS = ["KAPPA", "LAMP", "BACKGROUND",
                   "WAVE", "PINHOLE", "SLIT", "DARK", "FOIL"]


def get_spectral_binning(raw_dir) -> List[int]:
    """Gets the spectral binning according to the integration times used in the
    observation.

    Parameters
    ----------
    raw_dir : pathlib.Path
        The directory containing the raw observation files.

    Returns
    -------
    spectral_binning : list of int
        The spectral binning corresponding to the resolution of the
        observation.
    """
    readout = ReadoutFits(list(raw_dir.glob("*.fits"))[0])
    if readout.resolution == "high":
        resolution = f"{readout.resolution}_{readout.array_configuration}"
    else:
        resolution = readout.resolution
    return SPECTRAL_BINNING[resolution]


def set_script_arguments(mode: str) -> Tuple[str]:
    """Sets the arguments that are then passed to the
    "mat_autoPipeline" function.

    Parameters
    ----------
    raw_dir : pathlib.Path
        The directory containing the raw observation files.
    mode : str
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'.

    Returns
    -------
    lband_params : str
        The additional arguments passed to the `mat_autoPipeline` for the L-band. For the
        rest of the arguments see the `mat_autoPipeline`-script.
    nband_params : str
        The additional arguments passed to the `mat_autoPipeline` for the N-band. For the
        rest of the arguments see the `mat_autoPipeline`-script.

    See Also
    --------
    ..mat_tools.mat_autoPipeline.mat_autoPipeline
        A wrapper around esorex to execute the MATISSE-pipeline.
    """
    coh = "/corrFlux=TRUE/coherentAlgo=2/" if mode == "coherent" else ""
    return coh, f"{coh}/useOpdMod=TRUE/"


def prepare_reduction(raw_dir: Path,
                      calib_dir: Path,
                      product_dir: Path,
                      overwrite: Optional[bool]) -> None:
    """Prepares the reduction by removing removing old product files and
    sorting the raw files by associated calibrations and observations.

    Parameters
    ----------
    raw_dir : pathlib.Path
        The direcotry containing the raw observation files.
    calib_dir : pathlib.Path
        The directory containing the with the observation
        associated calibration files.
    product_dir : pathlib.Path
        The directory to contain the reduced files.
    overwrite : bool, optional
        If toggled, overwrites any files from previous reductions.
    """
    if not product_dir.exists():
        product_dir.mkdir(parents=True)
    if not calib_dir.exists():
        calib_dir.mkdir(parents=True)

    cprint("Checking files to move calibration files "
           "into 'calib_files' folders...", "g")
    for fits_file in raw_dir.glob("M.*"):
        shutil.move(fits_file, calib_dir / fits_file.name)

    # TODO: Make a good way to filter all the calibration files into
    # a folder directly
    # for fits_file in tqdm(list(raw_dir.glob("*.fits"))):
        # header = fits.getheader(fits_file, 0)
        # if not any(value in header["ESO DPR TYPE"]
                   # for value in ["OBJECT", "STD", "SKY"]):
            # shutil.move(fits_file, calib_dir / fits_file.name)


def cleanup_reduction(product_dir: Path,
                      mode: str, band: str,
                      do_data_quality_plot: bool,
                      maxIter: Optional[int] = None,
                      overwrite: Optional[bool] = None) -> None:
    """Moves the folders to their corresponding folders of structure
    "/mode/band" after the reduction has been finished and plots the
    (.fits)-files contained in them.

    Parameters
    ----------
    mode : str, optional
        The mode in which the reduction is to be executed. Either "coherent",
        "incoherent" or "both".
    band : str, optional
        The band in which the reduction is to be executed. Either "lband",
        "nband" or "both".
    do_data_quality_plot : bool, optional
        If toggled, plots the data quality of the incoherent L-band
        reduction.
    maxIter : int, optional
    overwrite : bool, optional
        If toggled, overwrites any files from previous reductions.
    """
    mode_and_band_dir = product_dir / mode / band
    if not mode_and_band_dir.exists():
        mode_and_band_dir.mkdir(parents=True)

    if maxIter == 1:
        iter_dir = "Iter1"
    else:
        if band == "nband":
            iter_dir = "Iter3"
        else:
            iter_dir = "Iter4"

    iter_dir = product_dir / iter_dir
    reduced_dirs = [folder for folder in iter_dir.iterdir() if folder.is_dir()]
    for reduced_folder in reduced_dirs:
        cprint(f"Copying folder '{reduced_folder.name}'...", "g")
        shutil.copytree(reduced_folder, mode_and_band_dir / reduced_folder.name,
                        dirs_exist_ok=True)

    cprint(f"Removing old iterations...", "g")
    for index in range(1, 6):
        if (product_dir / f"Iter{index}").exists():
            shutil.rmtree(product_dir / f"Iter{index}")

    reduced_dirs = [mode_and_band_dir / folder.name for folder in reduced_dirs]
    for reduced_folder in reduced_dirs:
        cprint(f"Plotting files of folder {reduced_folder.name}...", "g")
        for fits_file in get_fits_by_tag(reduced_folder, "RAW_INT"):
            plot_fits = Plotter(fits_file, save_path=reduced_folder)
            unwrap = True if "AQUARIUS" in str(fits_file) else False
            plot_fits.add_cphases(unwrap=unwrap).add_vis().add_vis2()
            plot_fits.plot(save=True, error=True)

        if do_data_quality_plot and mode == "incoherent" and band == "lband":
            plot_data_quality(
                    reduced_folder, reduced_folder / "data_quality")
    cprint(f"Finished reducing {band} in {mode}-mode", "lp")
    cprint(f"{'':-^50}", "lp")


def reduce_mode_and_band(raw_dir: Path, calib_dir: Path,
                         product_dir: Path, mode: bool,
                         band: str, ncores: int,
                         do_data_quality_plot: bool, overwrite: bool):
    """Reduces either the L- and/or the N-band data for either the 'coherent' and/or
    'incoherent' setting for a single iteration/epoch.

    Notes
    -----
    Removes the old (.sof)-files to ensure proper reduction and then creates
    folders in the 'res_dir'-directory. After this, it starts the reduction
    with the specified settings.

    Parameters
    ----------
    raw_dir : pathlib.Path
        The direcotry containing the raw observation files.
    calib_dir : pathlib.Path
        The directory containing the with the observation
        associated calibration files.
    product_dir : pathlib.Path
        The directory to contain the reduced files.
    mode : str
        The mode in which the reduction is to be executed. Either "coherent",
        "incoherent" or "both".
    band : str
        The band in which the reduction is to be executed. Either "lband",
        "nband" or "both".
    tpl_star : str
        The starting time of observations.
    ncores : int, optional
        The number of cores used.
    do_data_quality_plot : bool
        If toggled, plots the data quality of the incoherent L-band
        reduction.
    overwrite : bool
        If toggled, overwrites any files from previous reductions.
    """
    skip_L, skip_N = band == "nband", band == "lband"
    param_L, param_N = set_script_arguments(mode)
    spectral_binning = get_spectral_binning(raw_dir)

    # NOTE: here resol="" is required for the code not to skip the reduction
    reduction_kwargs = {"dirRaw": str(raw_dir), "dirResult": str(product_dir),
                        "dirCalib": str(calib_dir), "nbCore": ncores,
                        "paramL": param_L, "paramN": param_N, "resol": "",
                        "overwrite": int(overwrite), "maxIter": 1,
                        "skipL": int(skip_L), "skipN": int(skip_N),
                        "spectral_binning": spectral_binning}

    code = mat_autoPipeline(**reduction_kwargs)
    if code == -1:
        reduction_kwargs["maxIter"] = 5
        mat_autoPipeline(**reduction_kwargs)

    cleanup_reduction(product_dir, mode, band,
                      do_data_quality_plot,
                      reduction_kwargs["maxIter"], overwrite)


@print_execution_time
def reduction_pipeline(raw_dir: Path,
                       product_dir: Path,
                       mode: Optional[str] = "both",
                       band: Optional[str] = "both",
                       ncores: Optional[int] = 6,
                       do_data_quality_plot: Optional[bool] = True,
                       overwrite: Optional[bool] = False) -> None:
    """Runs the pipeline for the data reduction.

    Parameters
    ----------
    raw_dir : pathlib.Path
        The directory containing the raw observation files.
    product_dir : pathlib.Path
        The directory to contain the reduced files.
    mode : str, optional
        The mode in which the reduction is to be executed. Either "coherent",
        "incoherent" or "both".
    band : str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'.
    ncores : int, optional
        The number of cores used.
    do_data_quality_plot : bool, optional
        If toggled, plots the data quality of the incoherent L-band
        reduction.
    overwrite : bool, optional
        If toggled, overwrites any files from previous reductions.
    """
    raw_dir = Path(raw_dir).resolve()
    product_dir = Path(product_dir / "reduced").resolve()
    calib_dir = raw_dir / "calib_files"
    modes, bands = get_execution_modes(mode, band)

    prepare_reduction(raw_dir, calib_dir, product_dir, overwrite)
    for mode in modes:
        cprint(f"Processing the {mode} mode...", "lp")
        cprint(f"{'':-^50}", "lg")
        for band in bands:
            cprint(f"Processing the {band.title()}...", "lp")
            reduce_mode_and_band(raw_dir, calib_dir, product_dir,
                                 mode, band, ncores,
                                 do_data_quality_plot, overwrite)
    cprint(f"Finished reducing {', '.join(bands)} for {', '.join(modes)}-mode(s)", "lp")
