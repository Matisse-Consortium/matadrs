import shutil
import subprocess
from pathlib import Path
from typing import List, Optional
from collections import deque

from .fluxcal import fluxcal
from ..utils.plot import Plotter
from ..utils.tools import cprint, get_path_descriptor, check_if_target, get_fits_by_tag

__all__ = ["create_visbility_sof", "calibrate_files", "calibrate_folders", "calibrate"]

# TODO: Get correct calibrator for N- and L-band

DATA_DIR = Path(__file__).parent.parent.parent / "data"
DATABASE_DIR = DATA_DIR / "calibrator_databases"
DATABASES = ["vBoekelDatabase.fits", "calib_spec_db_v10.fits",
             "calib_spec_db_v10_supplement.fits", "calib_spec_db_supplement3.fits"]
LBAND_DATABASES = list(map(lambda x: DATABASE_DIR / x, DATABASES))
NBAND_DATABASES = LBAND_DATABASES[:]+[DATABASE_DIR / "vBoekelDatabase.fitsold"]

MODE_NAMES = {"coherent": "corrflux", "incoherent": "flux"}


# TODO: Make functionality that does not only calibrate willy nilly, but checks if the
# calibrator is there for N-band, L-band or LN-band, see Jozsef's files. Use the
# mat_target_list of newest edition for that (grep it via python api of google sheets)
# Make a function for this
# MAT_TARGET_LIST = DATA_DIR / "mat_target_list.xlsx"
def create_visbility_sof(reduced_dir: Path,
                         targets: List[Path], calibrators: List[Path]) -> Path:
    """Creates the (.sof)-file needed for the visibility calibration with 'mat_cal_oifits'
    and returns its path

    Parameters
    ----------
    reduced_dir: Path
        The direcotry containing the, by the pipeline, reduced data
    targets: List[Path]
        The 'TARGET_RAW_INT'-files
    calibrators: List[Path]
        The 'CALIB_RAW_INT'-files

    Returns
    -------
    sof_file: Path
        The Path to the (.sof)-file used for visibility calibration
    """
    sof_file = reduced_dir / "visbility_reduction.sof"
    with open(sof_file, "w+") as sof:
        for target in targets:
            sof.write(f"{target} TARGET_RAW_INT\n")
        sof.write("\n")
        for calibrator in calibrators:
            sof.write(f"{calibrator} CALIB_RAW_INT\n")
    return sof_file


# TODO: Implement visibility calibration like Jozsef, maybe?
# TODO: Make this better so it calibrates even if the calibrator or the science
# target is chopped but the other is not
def calibrate_files(reduced_dir: Path, target_dir: Path,
                    calibrator_dir: Path, mode: str) -> None:
    """The calibration for a target and a calibrator folder

    Parameters
    ----------
    reduced_dir: Path
        The path to multiple folders that need to be cross correlated
    target_dir: Path
        The directory of the science target
    calibrator_dir: Path
        The directory of the calibrator
    mode: str
        The mode of calibration. Either 'corrflux', 'flux' or 'both' depending
        if it is 'coherent' or 'incoherent' reduced data
    """
    # TODO: Make the following into a function
    cprint(f"Calibrating {target_dir.name} with {calibrator_dir.name}...", "p")
    targets = get_fits_by_tag(target_dir, "TARGET_RAW_INT")
    calibrators = get_fits_by_tag(calibrator_dir, "CALIB_RAW_INT")

    if not targets:
        cprint("No 'TARGET_RAW_INT*'-files found. SKIPPING!", "y")
        cprint(f"{'':-^50}", "lg")
        return
    if len(targets) != len(calibrators):
        cprint("#'TARGET_RAW_INT'-files != #'CALIB_RAW_INT'-files. SKIPPING!",
               "y")
        cprint(f"{'':-^50}", "lg")
        return

    output_dir = get_path_descriptor(reduced_dir, "TAR-CAL", targets[0], calibrators[0])
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    # TODO: Make the following into a function
    cprint(f"{'':-^50}", "lg")
    cprint("Calibrating fluxes...", "g")
    for index, (target, calibrator) in enumerate(zip(targets, calibrators), start=1):
        cprint(f"{'':-^50}", "lg")
        cprint(f"Processing {target.name} with {calibrator.name}...", "g")
        output_file = output_dir / f"TARGET_FLUX_CAL_INT_000{index}.fits"
        target, calibrator = str(target), str(calibrator)

        # TODO: Make this not as redundant, and the airmass correction implement as well?
        if "lband" in target:
            fluxcal(target, calibrator, str(output_file),
                    list(map(str, LBAND_DATABASES)), mode=MODE_NAMES[mode],
                    output_fig_dir=str(output_dir), do_airmass_correction=True)
        else:
            fluxcal(target, calibrator, str(output_file),
                    list(map(str, NBAND_DATABASES)), mode=MODE_NAMES[mode],
                    output_fig_dir=str(output_dir), do_airmass_correction=True)
        cprint("Plotting file...", "y")
        plot_fits = Plotter([output_file], save_path=output_dir)
        plot_fits.add_cphase().add_vis().plot(save=True)

    cprint(f"{'':-^50}", "lg")
    cprint("Calibrating visibilities...", "g")
    sof_file = create_visbility_sof(output_dir, targets, calibrators)
    subprocess.call([ESOREX_CMD, f"--output-dir={str(output_dir)}",
                     "mat_cal_oifits", str(sof_file)],
                    stdout=subprocess.DEVNULL)
    cprint("Plotting visibility files...", "g")
    for fits_file in get_fits_by_tag(output_dir, "TARGET_CAL_INT"):
        plot_fits = Plotter([fits_file], save_path=output_dir)
        plot_fits.add_cphase().add_vis().plot(save=True)

    # TODO: Find way to make this moving better than this -> Moves (.fits)-files
    # from this path that get created by 'mat_cal_oifits?'
    for fits_file in Path().cwd().glob("*.fits"):
        shutil.move(str(fits_file), str(output_dir / fits_file.name))
    shutil.move(str(Path().cwd() / "esorex.log"),
                str(output_dir / "mat_cal_oifits.log"))
    cprint(f"{'':-^50}", "lg")


def calibrate_folders(reduced_dir: Path, band: str, mode: str) -> None:
    """Takes two folders and calibrates their contents together

    Parameters
    ----------
    reduced_dir: Path
        The path to multiple folders that need to be cross correlated
    mode: str, optional
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'
    band: str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'
    """
    sub_dirs = sorted((reduced_dir / mode / band).glob("*.rb"))
    rotated_sub_directories = deque(sub_dirs.copy())
    rotated_sub_directories.rotate(1)

    for directory in sub_dirs:
        cprint(f"Calibration of {directory.name} in '{mode}'-mode'", "lp")
        cprint(f"{'':-^50}", "lg")
        if check_if_target(directory):
            for rotated_directory in rotated_sub_directories:
                calibrate_files(reduced_dir, directory, rotated_directory, mode)
        else:
            cprint("No 'TARGET_RAW_INT*'-files found. SKIPPED!", "y")
            cprint(f"{'':-^50}", "lg")
            continue


def calibrate(reduced_dir: Path, result_dir: Path,
              mode: Optional[str] = "both", band: Optional[str] = "both") -> None:
    """Does the full calibration for all of the "cal_dir" subdirectories

    Parameters
    ----------
    reduced_dir: Path
        The path to multiple folders that need to be cross correlated
    result_dir: Path
    mode: str, optional
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'
    band: str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'

    Notes
    -----
    The mode of calibration can either be 'coherent' or 'incoherent', which corresponds to
    either the calibration of the correlated fluxes or the total flux
    """
    if mode not in ["both", "coherent", "incoherent"]:
        raise IOError(f"No mode named '{mode}' exists!")
    if band not in ["both", "lband", "nband"]:
        raise IOError(f"No band named '{band}' exists!")
    modes = ["coherent", "incoherent"] if "both" else [mode]
    bands = ["lband", "nband"] if "both" else [band]

    for mode in modes:
        for band in bands:
            calibrate_folders(reduced_dir, band, mode)
    cprint("Calibration Done!", "lp")


    # TODO: Make this for every step of the reduction?
    # if remove_previous:
        # for directory in product_dir.glob("*"):
            # shutil.rmtree(directory)
        # cprint("Cleaning up old reduction...", "y")


if __name__ == "__main__":
    data_dir = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"
    stem_dir, target_dir = "matisse/GTO/hd142666/", "UTs/20220420"
    calibrate(data_dir, stem_dir, target_dir)
