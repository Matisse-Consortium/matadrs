import os

from glob import glob
from pathlib import Path
from warnings import warn
from collections import deque
from typing import Optional

from fluxcal import fluxcal

"""Smart wrapper for Jozsef Varga's script fluxcal.py"""

# Datapath to the directories to be calibrated, global variables
CAL_DATABASE_DIR = os.path.join(os.getcwd(), "calib_spec_databases")
CAL_DATABASE_FILES = ['vBoekelDatabase.fits', 'calib_spec_db_v10.fits',
                      'calib_spec_db_v10_supplement.fits']
CAL_DATABASE_PATHS = [os.path.join(CAL_DATABASE_DIR, i) for i in CAL_DATABASE_FILES]

# DATA_PATH = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"

# TODO: Add errorcodes to functions
# TODO: Also make the CAL-CAL calibration

def check_if_target(path: Path):
    """Checks if the given path contains TAR-files"""
    targets = glob(os.path.join(path, "TARGET_RAW_INT*"))
    return False if not targets else True

def get_path_descriptor(root_dir: Path, tar_dir: Path, cal_dir: Path) -> Path:
    """Assembles the names for the new directories that will contain the calibrated files
    and returns the "output_dir"

    Parameters
    ----------
    root_dir: Path
        The desired root directory
    tar_dir: Path
        The target directory
    cal_dir: Path
        The calibrator directory

    Returns
    -------
    output_dir: Path
        The path for the output directory
    """
    remote_working_dir = os.path.dirname(os.path.dirname(tar_dir))
    dir_name, time_sci, band = os.path.dirname(tar_dir).split('.')[:-1]
    dir_name = dir_name.split('/')[-1].replace("raw", "cal")
    time_cal = os.path.dirname(cal_dir).split('.')[-3]
    new_dir_name = '.'.join([root_dir, dir_name, time_sci, band, time_cal, "rb"])
    return os.path.join(remote_working_dir, "calib", new_dir_name)

def single_calibration(tar_dir: Path, cal_dir: Path, mode: str) -> None:
    """The calibration for a target and a calibrator folder

    Parameters
    ----------
    tar_dir: Path
        A specific directory to be the target for calibration
    cal_dir: Path
        A specific directory to be the calibrator for calibration
    mode: str
        The mode of calibration. Either "corrflux", "flux" or "both" depending
        if it is "coherent" or "incoherent" reduced data

    See Also
    --------
    get_path_descriptor()
    """
    print(f"Calibrating {os.path.basename(tar_dir)} with "\
          f"{os.path.basename(cal_dir)}")
    targets = glob(os.path.join(tar_dir, "TARGET_RAW_INT*"))
    if not targets:
        warn("No 'TARGET_RAW_INT*'-files found. SKIPPED!")
        print("------------------------------------------------------------")
        return

    targets.sort(key=lambda x: x[-8:])
    calibrators = glob(os.path.join(cal_dir, "CALIB_RAW_INT*"))
    calibrators.sort(key=lambda x: x[-8:])

    if len(targets) != len(calibrators):
        warn("#'TARGET_RAW_INT'-files != #'CALIB_RAW_INT'-files. SKIPPING!")
        print("------------------------------------------------------------")
        return

    output_dir = get_path_descriptor("TAR-CAL", targets[0], calibrators[0])
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, target in enumerate(targets):
        print("------------------------------------------------------------")
        print(f"Processing {os.path.basename(target)} with "\
              f"{os.path.basename(calibrators[i])}")
        output_file = os.path.join(output_dir, f"TARGET_CAL_INT_000{i+1}.fits")

        fluxcal(target, calibrators[i], output_file,\
                CAL_DATABASE_PATHS, mode=mode, output_fig_dir=output_dir)
    print("------------------------------------------------------------")
    print("Done!")
    print("------------------------------------------------------------")

def do_calibration(root_dir: Path,
                   mode: Optional[str] = "corrflux") -> None:
    """Takes two folders and calibrates their contents together

    Parameters
    ----------
    root_dir: Path
        The path to multiple folders that need to be cross correlated. Will be
        skipped if folders for targets and calibrators are specified
    mode: str, optional
        The mode of calibration. Either "corrflux", "flux" or "both" depending
        if it is "coherent" or "incoherent". Default mode is "corrflux"
    """
    subdirs = glob(os.path.join(root_dir, "*.rb"))
    subdirs.sort(key=lambda x: x.split(".")[~2])
    subdirs_rotated = deque(subdirs.copy())
    subdirs_rotated.rotate(1)

    for path in subdirs:
        if check_if_target(path):
            for path_rotated in subdirs_rotated:
                paths = (path, path_rotated)
                single_calibration(*paths, mode=mode)
        else:
            warn("No 'TARGET_RAW_INT*'-files found. SKIPPED!")
            print("------------------------------------------------------------")
            continue

def calibration_pipeline(reduced_dir: Path, both: Optional[bool] = True,
                         lband: Optional[bool] = False) -> None:
    """Does the full calibration for all of the "cal_dir" subdirectories

    Parameters
    ----------
    reducd_dir: Path
        The folder to be calibrated. Subdirectories will automatically be calibrated
    both: bool, optional
        If both bands are to be calibrated, this has to be false for the "lband" option
        to be considered
    lband: bool, optional
        If "both=False" and this is "True"", then lband will be calibrated, if
        "both=False" and this is "False", then nband will be calibrated
    """
    modes, bands = {"coherent": "corrflux", "incoherent": "flux"},\
            ["lband", "nband"]

    if both:
        for i, o in modes.items():
            path = os.path.join(reduced_dir, i)
            for j in bands:
                temp_path = os.path.join(path, j)
                print(f"Calibration of {temp_path}")
                print(f"with mode={o}")
                print("------------------------------------------------------------")
                do_calibration(temp_path, mode=o)
    else:
        for i, o in modes.items():
            band = "lband" if lband else "nband"
            path = os.path.join(reduced_dir, i, band)
            print(f"Calibration of {path}")
            print(f"with mode={o}")
            print("------------------------------------------------------------")
            do_calibration(path, mode=o)


if __name__ == "__main__":
    specific_dir = "/Users/scheuck/Documents/data/20190514"
    # specific_path = "GTO/hd142666/PRODUCT/UTs/20190514"
    calibration_pipeline(specific_dir, both=False)
    # calibration_pipeline(os.path.join(DATA_PATH, specific_path), both=True)

