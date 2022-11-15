import os
import time
import datetime

from warnings import warn
from pathlib import Path
from typing import List, Optional

import mat_autoPipeline as mp

# TODO: Look up high spectral binning and make savefile somehow show all
# High spectral binning is 7, 49


def set_script_arguments(do_corr_flux: bool, array: str,
                         lband: Optional[bool] = False,
                         spectral_binning: Optional[List] = [5, 7]) -> str:
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
    binning_L, binning_N = spectral_binning
    tel = 3 if array == "ATs" else 0
    opd_mod = "useOpdMod=FALSE/" if lband else "useOpdMod=TRUE/"
    flux = f"corrFlux=TRUE/{opd_mod}coherentAlgo=2/" if do_corr_flux else ""
    paramL_lst = f"/spectralBinning={binning_L}/{flux}compensate='pb,rb,nl,if,bp,od'"
    paramN_lst = f"/replaceTel={tel}/{flux}spectralBinning={binning_N}"
    return paramL_lst, paramN_lst


def single_reduction(raw_dir: str, calib_dir: str, res_dir: str,
                     array: str, mode: bool, band: str) -> None:
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
    start_time = time.time()
    path_lst = ["coherent" if mode else "incoherent", "lband" if band else "nband"]
    path = "/".join(path_lst)
    sub_dir = os.path.join(res_dir, path)
    paramL, paramN = set_script_arguments(mode, array, band)
    skipL, skipN = int(not band), int(band)

    # NOTE: Removes the old '.sof'-files
    try:
        os.system(f"rm {os.path.join(res_dir, 'Iter1/*.sof*')}")
        os.system(f"rm -r {os.path.join(res_dir, 'Iter1/*.rb')}")
        os.system(f"rm -r {os.path.join(sub_dir, '*.rb')}")
        print("Old (.sof)-files have been deleted!")
    # TODO: Make logger here
    except Exception:
        warn(f"Removing of (.sof)- and (.rb)-files from {sub_dir} has failed!")

    if not os.path.exists(sub_dir):
        os.makedirs(sub_dir)

    mp.mat_autoPipeline(dirRaw=raw_dir, dirResult=res_dir,
                        dirCalib=calib_dir,
                        nbCore=6, resol='',
                        paramL=paramL, paramN=paramN,
                        overwrite=0, maxIter=1,
                        skipL=skipL, skipN=skipN)

    try:
        os.system(f"mv -f {os.path.join(res_dir, 'Iter1/*.rb')} {sub_dir}")
    # TODO: Make logger here
    except Exception:
        print("Moving of files to {sub_dir} failed!")

    # Takes the time at end of execution
    print("---------------------------------------------------------------------")
    print(f"Executed the {path_lst[0]} {path_lst[1]} reduction in"
          f" {datetime.timedelta(seconds=(time.time()-start_time))} hh:mm:ss")
    print("---------------------------------------------------------------------")


def reduction_pipeline(raw_dir: Path, calib_dir: Path, res_dir: Path,
                       array: str, both: Optional[bool] = False,
                       lband: Optional[bool] = False):
    """Runs the pipeline for the data reduction

    Parameters
    ----------
    raw_dir: Path
        The path containing the raw observation files
    calib_dir: Path
        The path containing the calibration files
    res_dir: Path
        The path to contain to reduced data
    array: str
        The array configuration that was used for the observation
    both: bool, optional
        If both bands are to be reduced, this has to be false for the "lband" option to be
        considered
    lband: bool, optional
        If "both=False" and this is "True"", then lband will be calibrated, if
        "both=False" and this is "False", then nband will be calibrated

    See Also
    --------
    single_reduction()
    """
    overall_start_time = time.time()
    if both:
        band_bools = [True, False]
    else:
        band_bools = [lband]

    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    for mode_bools in [True, False]:
        mode = "coherent" if mode_bools else "incoherent"
        print(f"Processing {mode} reduction")
        print("---------------------------------------------------------------------")
        for band in band_bools:
            single_reduction(raw_dir, calib_dir, res_dir, array,\
                             mode=mode_bools, band=band)

    print(f"Executed the overall reduction in"
          f" {datetime.timedelta(seconds=(time.time()-overall_start_time))} hh:mm:ss")


if __name__ == "__main__":
    data_path = "/data/beegfs/astro-storage/groups/matisse/scheuck/data"
    stem_dir, target_dir = "matisse/GTO/hd163296/", "ATs/20190323"
    raw_dir = calib_dir = os.path.join(data_path, stem_dir, "RAW", target_dir)
    res_dir = os.path.join(data_path, stem_dir, "PRODUCTS", target_dir)
    reduction_pipeline(raw_dir, calib_dir, res_dir, "ATs", both=True)
