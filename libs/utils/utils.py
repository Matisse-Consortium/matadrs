import os

from glob import glob
from pathlib import Path
from shutil import copyfile
from typing import Optional, List

from astropy.io import fits


def cprint(message: str, c: Optional[str] = None) -> None:
    """Prints with color"""
    color_dict = {"r": ["\033[91m", "\033[00m"], "g": ["\033[92m", "\033[00m"],
                  "y": ["\033[93m", "\033[00m"], "lp": ["\033[94m", "\033[00m"],
                  "p": ["\033[95m", "\033[00m"], "cy": ["\033[96m", "\033[00m"],
                  "lg": ["\033[97m", "\033[00m"]}

    if c:
        colored_string = color_dict[c]
        colored_string.insert(1, message)
        print("".join(colored_string))
    else:
        print(message)

def split_fits(folder: Path, tag: str):
    """Searches a folder for a tag and returns the non-chopped
    and chopped (.fits)-files"""
    unchopped_fits = get_fits_by_tag(folder, tag)
    if len(unchopped_fits) == 6:
        return unchopped_fits[:4], unchopped_fits[4:]
    return unchopped_fits, None


def get_fits_by_tag(folder: Path, tag: str):
    """Searches a folder for a tag and returns the (.fits)-files matching it"""
    return sorted(folder.glob(f"*{tag}*.fits"), key=lambda x: x.name[-8:])


def check_if_target(target_dir: Path) -> bool:
    """Checks if the given path contains TAR-files

    Parameters
    ----------
    target_dir: Path
        The directory that is to be checked for 'TARGET_RAW_INT' files

    Returns
    -------
    contains_target: bool
    """
    return True if glob(os.path.join(target_dir, "TARGET_RAW_INT*")) else False


def get_path_descriptor(root_dir: Path, descriptor: Path,
                        tar_dir: Path, cal_dir: Path) -> Path:
    """Assembles the names for the new directories that will contain the
    calibrated files and returns the 'output_dir'

    Parameters
    ----------
    root_dir: Path
        The root directory of the PRODUCT
    descriptor: Path
        The desired name 'TAR-CAL' directory
    tar_dir: Path
        The target directory
    cal_dir: Path
        The calibrator directory

    Returns
    -------
    output_dir: Path
        The path for the output directory
    """
    mode_and_band = str(tar_dir.parents[1]).split("/")[-2:]
    dir_name, time_stamp_sci, detector = str(tar_dir.parent).split(".")[:-1]
    dir_name = dir_name.split('/')[-1].replace("raw", "cal")
    time_stamp_cal = str(cal_dir.parent).split('.')[-3]
    new_dir_name = '.'.join([descriptor, dir_name,
                             time_stamp_sci, detector, time_stamp_cal, "rb"])
    return root_dir / "calib" / mode_and_band[0] / new_dir_name


if __name__ == "__main__":
    print(cprint("Hello", "lp"))
