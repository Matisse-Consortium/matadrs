import os

from glob import glob
from pathlib import Path
from typing import Optional


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
    mode_and_band = os.path.dirname(os.path.dirname(tar_dir)).split("/")[-2:]
    dir_name, time_stamp_sci, band = os.path.dirname(tar_dir).split('.')[:-1]
    dir_name = dir_name.split('/')[-1].replace("raw", "cal")
    time_cal = os.path.dirname(cal_dir).split('.')[-3]
    new_dir_name = '.'.join([descriptor, dir_name,
                             time_stamp_sci, band, time_cal, "rb"])
    return os.path.join(root_dir, "calib", mode_and_band[0], new_dir_name)


if __name__ == "__main__":
    print(cprint("Hello", "lp"))
