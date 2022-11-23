import os
import time
import datetime
import shutil

from glob import glob
from pathlib import Path
from typing import List, Optional

import mat_autoPipeline as mp

from utils import cprint

SPECTRAL_BINNING = {"LOW": [5, 7], "HIGH": [7, 49]}

def set_script_arguments(do_corr_flux: bool, array: str,
                         resolution: Optional[str] = "LOW") -> str:
    """Sets the arguments that are then passed to the 'mat_autoPipeline.py'
    script

    Parameters
    ----------
    do_corr_flux: bool
        This specifies if the flux is to be reduced or not
    array: str
        The array configuration that was used for the observation
    spectral_binning: List, optional
        The spectral "binning" to be selected

    Returns
    -------
    str
        A string that contains the arguments, which are passed to the MATISSE-pipline
    """
    bin_L, bin_N = SPECTRAL_BINNING[resolution]
    # NOTE: Jozsef uses 3 here, but Jacob 2? What is the meaning of this -> Read up on it
    # -> Already asked, awaiting response
    tel = "/replaceTel=3" if array == "ATs" else ""
    coh_L  = f"corrFlux=TRUE/coherentAlgo=2/" if do_corr_flux else ""
    coh_N = f"corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2/"\
            if do_corr_flux else ""
    paramL_lst = f"/{coh_L}compensate=[pb,rb,nl,if,bp,od]/spectralBinning={bin_L}"
    paramN_lst = f"/{coh_N}compensate=[pb,rb,nl,if,bp,od]/spectralBinning={bin_N}{tel}"
    return paramL_lst, paramN_lst


def single_reduction(raw_dir: Path, res_dir: Path,
                     array: str, mode: bool, band: str,
                     resolution: Optional[str] = "LOW") -> None:
    """Reduces either the lband or the nband data for either the "coherent" or
    "incoherent" setting for a single iteration/epoch.

    Removes the old (.sof)-files to ensure proper reduction and then creates the needed
    folders in the "res_dir"-directory. After this, it starts the reduction with the
    specified settings

    Parameters
    ----------
    raw_dir: Path
        The path containing the raw observation files
    res_dir: Path
        The path to contain to reduced data
    array: str
        The array configuration that was used for the observation. Either "AT" or "UT"
    mode: bool
        The mode in which the reduction is to be done, either "incoherent" if "False" or
        "coherent" if "True"
    band: str
        The band for which the reduction is to be done, either "lband" or "nband"

    See Also
    --------
    set_script_arguments()
    """
    # TODO: Replace this time with process finish time
    start_time = time.time()
    mode_and_band_dir = os.path.join(res_dir, mode, band)
    param_L, param_N = set_script_arguments(mode, array, resolution)
    skip_L = True if band == "nband" else False
    skip_N = not skip_L

    try:
        shutil.rmtree(os.path.join(res_dir, "Iter1"))
        cprint("Old 'Iter1'-folder has been deleted!", "g")
    # TODO: Make logger here
    except Exception:
        cprint(f"Removing of 'Iter1' folder from {mode_and_band_dir} has failed!",
               "y")

    if not os.path.exists(mode_and_band_dir):
        os.makedirs(mode_and_band_dir)

    # mp.mat_autoPipeline(dirRaw=raw_dir, dirResult=res_dir, dirCalib=raw_dir,
                        # nbCore=6, resol='', paramL=param_L, paramN=param_N,
                        # overwrite=0, maxIter=1, skipL=skip_L, skipN=skip_N)

    try:
        rb_folders = glob(os.path.join(res_dir, "Iter1/*.rb"))
        for folder in rb_folders:
            if os.path.exists(os.path.join(mode_and_band_dir,
                                           os.path.basename(folder))):
                shutil.rmtree(os.path.join(mode_and_band_dir,
                                           os.path.basename(folder)))
            shutil.move(folder, mode_and_band_dir)
        cprint("Folders have sucessfully been moved to their directories", "g")
    # TODO: Make logger here
    except Exception:
        cprint("Moving of files to {mode_and_band_dir} failed!", "y")

    print("---------------------------------------------------------------------")
    cprint(f"Executed the {mode} reduction for the {band} in"
          f" {datetime.timedelta(seconds=(time.time()-start_time))} hh:mm:ss")
    cprint("---------------------------------------------------------------------",
          "lg")


def reduction_pipeline(root_dir: Path, stem_dir: Path,
                       target_dir: Path, array: str,
                       resolution: Optional[str] = "LOW"):
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
    # TODO: Replace this time with process finish time
    overall_start_time = time.time()
    raw_dir = os.path.join(root_dir, stem_dir, "RAW", target_dir)
    res_dir = os.path.join(root_dir, stem_dir, "PRODUCTS", target_dir)

    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    for mode in ["coherent", "incoherent"]:
        cprint(f"Processing {mode} reduction", "lp")
        cprint("---------------------------------------------------------------------",
              "lg")
        for band in ["lband", "nband"]:
            single_reduction(raw_dir, res_dir, array,
                             mode=mode, band=band)

    cprint(f"Executed the overall reduction in"
           f" {datetime.timedelta(seconds=(time.time()-overall_start_time))} hh:mm:ss",
          "lp")


if __name__ == "__main__":
    data_dir = "/data/beegfs/astro-storage/groups/matisse/scheuck/data"
    stem_dir, target_dir = "matisse/GTO/hd163296/", "ATs/20190323"
    reduction_pipeline(data_dir, stem_dir, target_dir, "ATs")
