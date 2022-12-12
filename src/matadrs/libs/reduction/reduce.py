import time
import datetime
import shutil

from pathlib import Path
from typing import Optional

from mat_tools import mat_autoPipeline as mp

from .utils import cprint


SPECTRAL_BINNING = {"low": [5, 7], "high_ut": [5, 38], "high_at": [5, 98]}


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
    # NOTE: Jozsef uses 3 here, but Jacob 2? What is the meaning of this -> Read up on it
    # -> Already asked, awaiting response
    compensate = "/compensate=[pb,rb,nl,if,bp,od]"
    tel = "/replaceTel=3" if array == "ATs" else "/replaceTel=0"
    coh_L  = f"/corrFlux=TRUE/useOpdMod=FALSE/coherentAlgo=2"\
            if corr_flux else ""
    coh_N = f"/corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2" if corr_flux else ""
    paramL_lst = f"{coh_L}{compensate}/spectralBinning={bin_L}"
    paramN_lst = f"{tel}{coh_N}/spectralBinning={bin_N}"
    return paramL_lst, paramN_lst


def reduce_mode_and_band(raw_dir: Path, calib_dir: Path, res_dir: Path,
                         array: str, mode: bool, band: str,
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

    See Also
    --------
    set_script_arguments()
    """
    # TODO: Replace this time with process finish time
    start_time = time.time()
    mode_and_band_dir = res_dir / mode / band
    param_L, param_N = set_script_arguments(mode, array, resolution)
    skip_L = True if band == "nband" else False
    skip_N = not skip_L

    if not mode_and_band_dir.exists():
        mode_and_band_dir.mkdir()

    mp.mat_autoPipeline(dirRaw=raw_dir, dirResult=res_dir, dirCalib=calib_dir,
                        nbCore=6, resol='', paramL=param_L, paramN=param_N,
                        overwrite=0, maxIter=1, skipL=skip_L, skipN=skip_N)

    try:
        rb_folders = res_dir.glob("Iter1/*.rb")
        for folder in rb_folders:
            if (mode_and_band_dir / folder.name).exists():
                shutil.rmtree(mode_and_band_dir / folder.name)
            shutil.move(folder, mode_and_band_dir)

        if rb_folders:
            cprint("Folders have sucessfully been moved to their directories",
                   "g")
    # TODO: Make logger here
    except Exception:
        cprint("Moving of files to {mode_and_band_dir} failed!", "y")

    print("---------------------------------------------------------------------")
    cprint(f"Executed the {mode} reduction for the {band} in"
          f" {datetime.timedelta(seconds=(time.time()-start_time))} hh:mm:ss")
    cprint("---------------------------------------------------------------------",
          "lg")


def reduce(root_dir: Path, stem_dir: Path,
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
    # TODO: Replace this time with process finish time
    overall_start_time = time.time()
    raw_dir = Path(root_dir, stem_dir, "raw", target_dir)

    # TODO: Change this to proper search for calibration_files
    calib_dir = raw_dir
    res_dir = Path(root_dir, stem_dir, "products", target_dir)

    if not res_dir.exists():
        res_dir.mkdir()

    # TODO: Add in the option to not remove old reduction and make new one take an
    # addional tag after its name
    try:
        folders = res_dir.glob("*")
        for folder in folders:
            if folder.exists():
                shutil.rmtree(folder)
        cprint("Cleaned up old reduction!", "y")

    # TODO: Make logger here
    except Exception:
        cprint("Cleaning up failed!", "y")


    for mode in ["coherent", "incoherent"]:
        cprint(f"Processing {mode} reduction", "lp")
        cprint("---------------------------------------------------------------------",
              "lg")
        for band in ["lband", "nband"]:
            reduce_mode_and_band(raw_dir, calib_dir, res_dir, array, mode=mode, band=band)
            breakpoint()

    cprint(f"Executed the overall reduction in"
           f" {datetime.timedelta(seconds=(time.time()-overall_start_time))} hh:mm:ss",
          "lp")


if __name__ == "__main__":
    data_dir = "../../../../data"
    stem_dir, target_dir = "/hd163296/", "ATs/20190323"
    reduce(data_dir, stem_dir, target_dir, "ATs")

