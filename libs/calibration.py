import os
import warnings
import matplotlib.pyplot as plt

from glob import glob
from pathlib import Path
from collections import deque
from typing import Optional

from plot import Plotter
from utils import get_path_descriptor, check_if_target
from fluxcal import fluxcal

# Datapath to the directories to be calibrated, global variables
DATABASE_DIR = os.path.join(os.path.dirname(__file__), "calib_spec_databases")
DATABASE_FILES = ['vBoekelDatabase.fits', 'calib_spec_db_v10.fits',
                  'calib_spec_db_v10_supplement.fits']
DATABASE_PATHS = [os.path.join(DATABASE_DIR, databases) for databases in DATABASE_FILES]

# TODO: Also make the CAL-CAL calibration


def single_calibration(root_dir: Path, tar_dir: Path,
                       cal_dir: Path, mode_name: str) -> None:
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
        warnings.warn("No 'TARGET_RAW_INT*'-files found. SKIPPED!")
        print("------------------------------------------------------------")
        return

    targets.sort(key=lambda x: x[-8:])
    calibrators = glob(os.path.join(cal_dir, "CALIB_RAW_INT*"))
    calibrators.sort(key=lambda x: x[-8:])

    if len(targets) != len(calibrators):
        warnings.warn("#'TARGET_RAW_INT'-files != #'CALIB_RAW_INT'-files. SKIPPING!")
        print("------------------------------------------------------------")
        return

    output_dir = get_path_descriptor(root_dir, "TAR-CAL",
                                     targets[0], calibrators[0])
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # TODO: Fix the numbering of the (.fits)-files
    for index, (target, calibrator) in enumerate(zip(targets, calibrators)):
        print("------------------------------------------------------------")
        print(f"Processing {os.path.basename(target)} with "\
              f"{os.path.basename(calibrator)}")
        output_file = os.path.join(output_dir, f"TARGET_CAL_INT_000{index+1}.fits")

        fluxcal(target, calibrator, output_file,\
                DATABASE_PATHS, mode=mode_name,
                output_fig_dir=output_dir)

    print("------------------------------------------------------------")
    print("Creating plots...")
    fits_files = glob(os.path.join(output_dir, "*.fits"))
    for fits_file in fits_files:
        plot_fits = Plotter([fits_file], save_path=output_dir)
        plot_fits.add_cphases().add_corr_flux().plot(save=True)
    print("Plots created!")
    print("------------------------------------------------------------")
    print("Done!")
    print("------------------------------------------------------------")


def do_calibration(root_dir: Path, band_dir: Path,
                   mode_name: Optional[str] = "corrflux") -> None:
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
    sub_dirs = glob(os.path.join(os.path.join(root_dir, band_dir), "*.rb"))
    sub_dirs.sort(key=lambda x: x.split(".")[~2])
    sub_dirs_rotated = deque(sub_dirs.copy())
    sub_dirs_rotated.rotate(1)

    for directory in sub_dirs:
        if check_if_target(directory):
            for dir_rotated in sub_dirs_rotated:
                single_calibration(root_dir, directory,
                                   dir_rotated, mode_name=mode_name)
        else:
            warnings.warn("No 'TARGET_RAW_INT*'-files found. SKIPPED!")
            print("------------------------------------------------------------")
            continue


def calibration_pipeline(data_dir: Path, stem_dir: Path, target_dir: Path):
    """Does the full calibration for all of the "cal_dir" subdirectories

    Parameters
    ----------
    data_dir: Path
    stem_dir: Path
    target_dir: Path
    """
    root_dir = os.path.join(data_dir, stem_dir, "PRODUCTS", target_dir)
    modes, bands = {"coherent": "corrflux", "incoherent": "flux"},\
            ["lband", "nband"]

    for mode, mode_name in modes.items():
        for band in bands:
            band_dir = os.path.join(mode, band)
            print(f"Calibration of {band_dir}")
            print(f"with mode_name={mode_name}")
            print("------------------------------------------------------------")
            do_calibration(root_dir, band_dir, mode_name=mode_name)


if __name__ == "__main__":
    data_dir = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"
    stem_dir, target_dir = "matisse/GTO/hd163296/", "ATs/20190323"
    calibration_pipeline(data_dir, stem_dir, target_dir)

