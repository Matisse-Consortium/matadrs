import time
from datetime import timedelta
from functools import wraps
from pathlib import Path
from typing import Callable, Tuple, List, Optional

__all__ = ["cprint", "print_execution_time", "get_execution_modes",
           "split_fits", "get_fits_by_tag", "check_if_target", "get_path_descriptor"]


def cprint(message: str, color: Optional[str] = None) -> None:
    """Makes use of ascii-codes to print messages in color
    Parameters
    ----------
    message: str
        The message to be printed
    color: str
        The name of the color to be used
    """
    color_dict = {"r": ["\033[91m", "\033[00m"], "g": ["\033[92m", "\033[00m"],
                  "y": ["\033[93m", "\033[00m"], "lp": ["\033[94m", "\033[00m"],
                  "p": ["\033[95m", "\033[00m"], "cy": ["\033[96m", "\033[00m"],
                  "lg": ["\033[97m", "\033[00m"]}

    if color is not None:
        color_code = color_dict[color]
        print(color_code[0] + message + color_code[1])
    else:
        print(message)


def print_execution_time(func: Callable):
    """Prints the execution time of the input function"""
    @wraps(func)
    def inner(*args, **kwargs):
        overall_start_time = time.perf_counter()
        result = func(*args, **kwargs)
        execution_time = time.perf_counter()-overall_start_time
        cprint(f"Executed in {timedelta(seconds=execution_time)}"
               " hh:mm:ss", "lg")
        return result
    return inner


def get_execution_modes(mode: Optional[str] = None,
                        band: Optional[str] = None) -> Tuple[List[str], List[str]]:
    """Determines the mode- and band configurations used by the users input. Returns
    either one or two lists depending on the input

    Parameters
    ----------
    mode: str, optional
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'
    band: str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'

    Returns
    -------
    modes: List[str]
        A list of the modes to be enhanced
    bands: List[str]
        A list of the bands to be enhanced
    """
    modes, bands = [], []
    if mode is not None:
        if mode not in ["both", "coherent", "incoherent"]:
            raise IOError(f"No mode named '{mode}' exists!")
        modes = ["coherent", "incoherent"] if "both" else [mode]
    if band is not None:
        if band not in ["both", "lband", "nband"]:
            raise IOError(f"No band named '{band}' exists!")
        bands = ["lband", "nband"] if "both" else [band]
    return modes, bands


def split_fits(directory: Path, tag: str) -> Tuple[List[Path], Optional[List[Path]]]:
    """Searches a folder for a tag and if files are found it returns the non-chopped
    and chopped (.fits)-files. If there are only non-chopped (.fits)-files it will return
    'None' for the chopped-files

    Parameters
    ----------
    directory: Path
        The directory to be searched in
    tag: str
        The tag that is contained in the file names

    Returns
    ----------
    unchopped_fits: List[Path]
        A list of Paths that are the chopped (.fits)-files
    chopped_fits: List[Path] | None
        A list of Paths that are the unchopped (.fits)-files
    """
    unchopped_fits = get_fits_by_tag(directory, tag)
    if len(unchopped_fits) == 6:
        return unchopped_fits[:4], unchopped_fits[4:]
    return unchopped_fits, None


def get_fits_by_tag(directory: Path, tag: str) -> List[Path]:
    """Searches a folder for a tag and returns the (.fits)-files matching it

    Parameters
    ----------
    directory: Path
        The directory to be searched in
    tag: str
        The tag that is contained in the file names

    Returns
    ----------
    files: List[Path]
        A list of files that contain the tag in their names
    """
    return sorted(directory.glob(f"*{tag}*.fits"), key=lambda x: x.name[-8:])


def check_if_target(target_dir: Path) -> bool:
    """Checks if the given directory contains 'TARGET_RAW_INT'-files

    Parameters
    ----------
    target_dir: Path
        The directory that is to be checked for 'TARGET_RAW_INT' files

    Returns
    -------
    contains_target: bool
    """
    return True if target_dir.glob("TARGET_RAW_INT*") else False


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
