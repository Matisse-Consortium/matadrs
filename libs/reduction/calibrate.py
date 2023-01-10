import shutil
import subprocess
from pathlib import Path
from typing import List, Optional
from collections import deque, namedtuple

# TODO: Find way to make this into a complete module -> More pythonic!
from plot import Plotter
from readout import ReadoutFits
from fluxcal import fluxcal
from calib_BCD2 import calib_BCD
from utils import get_path_descriptor, check_if_target, cprint


# TODO: Also make the CAL-CAL calibration??
DATA_DIR = Path(__file__).parent.parent.parent / "data"
DATABASE_DIR = DATA_DIR / "calibrator_databases"
LBAND_DATABASES = [DATABASE_DIR / database\
        for database in ["vBoekelDatabase.fits", "calib_spec_db_v10.fits",
                        "calib_spec_db_v10_supplement.fits",
                         "calib_spec_db_supplement3.fits"]]
NBAND_DATABASES = LBAND_DATABASES[:]+[DATABASE_DIR / "vBoekelDatabase.fitsold"]

ESOREX_CMD = "/data/beegfs/astro-storage/groups/matisse/isbell/esorex_installation/bin/esorex"


# TODO: Make functionality that does not only calibrate willy nilly, but checks if the
# calibrator is there for N-band, L-band or LN-band, see Jozsef's files. Use the
# mat_target_list of newest edition for that (grep it via python api of google sheets)
# Make a function for this
# MAT_TARGET_LIST = DATA_DIR / "mat_target_list.xlsx"
def create_visbility_sof(raw_dir: Path,
                         targets: List[Path], calibrators: List[Path]) -> Path:
    """Creates the (.sof)-file needed for the visibility calibration with mat_cal_oifits"""
    sof_file = raw_dir / "visbility_reduction.sof"
    with open(sof_file, "w+") as sof:
        for target in targets:
            sof.write(f"{target} TARGET_RAW_INT\n")
        sof.write("\n")
        for calibrator in calibrators:
            sof.write(f"{calibrator} CALIB_RAW_INT\n")
    return sof_file


# TODO: Implement visibility calibration like Jozsef, maybe?
def calibrate_fits_files(root_dir: Path, tar_dir: Path,
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
    # TODO: Make the following into a function
    cprint(f"Calibrating {tar_dir.name} with {cal_dir.name}...", "p")
    targets = sorted(tar_dir.glob("TARGET_RAW_INT*"),
                     key=lambda x: x.name[-8:])

    if not targets:
        cprint("No 'TARGET_RAW_INT*'-files found. SKIPPED!", "y")
        cprint(f"{'':-^50}", "lg")
        return

    calibrators = sorted(cal_dir.glob("CALIB_RAW_INT*"),
                         key=lambda x: x.name[-8:])

    # TODO: Make this better so it calibrates even if the calibrator or the science
    # target is chopped but the other is not
    if len(targets) != len(calibrators):
        cprint("#'TARGET_RAW_INT'-files != #'CALIB_RAW_INT'-files. SKIPPING!", "y")
        cprint(f"{'':-^50}", "lg")
        return

    output_dir = get_path_descriptor(root_dir, "TAR-CAL", targets[0], calibrators[0])
    mode = "incoherent" if "incoherent" in str(output_dir) else "coherent"
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    # TODO: Make the following into a function
    cprint(f"{'':-^50}", "lg")
    cprint("Calibrating fluxes...", "g")
    for index, (target, calibrator) in enumerate(zip(targets, calibrators), start=1):
        cprint(f"{'':-^50}", "lg")
        cprint(f"Processing {target.name} with {calibrator.name}...", "g")
        output_file = str(output_dir / f"TARGET_FLUXCAL_INT_000{index}.fits")
        target, calibrator = str(target), str(calibrator)

        # TODO: Make this not as redundant, and the airmass correction implement as well?
        if "lband" in target:
            fluxcal(target, calibrator, output_file,
                    list(map(str, LBAND_DATABASES)), mode=mode_name,
                    output_fig_dir=str(output_dir), do_airmass_correction=True)
        else:
            fluxcal(target, calibrator, output_file,
                    list(map(str, NBAND_DATABASES)), mode=mode_name,
                    output_fig_dir=str(output_dir), do_airmass_correction=True)

    cprint(f"{'':-^50}", "lg")
    cprint("Calibrating visibilities...", "g")
    sof_file = create_visbility_sof(output_dir, targets, calibrators)
    subprocess.call([ESOREX_CMD, f"--output-dir={str(output_dir)}",
                     "mat_cal_oifits", str(sof_file)],
                    stdout=subprocess.DEVNULL)

    # TODO: Find way to make this moving better than this -> Moves (.fits)-files
    # from this path that get created by 'mat_cal_oifits?'
    for fits_file in Path().cwd().glob("*.fits"):
        shutil.move(str(fits_file), str(output_dir / fits_file.name))
    shutil.move(str(Path().cwd() / "esorex.log"),
                str(output_dir / "mat_cal_oifits.log"))

    # TODO: Fix this at some point
    cprint(f"{'':-^50}", "lg")
    cprint("Creating plots...", "g")
    for fits_file in output_dir.glob("*.fits"):
        plot_fits = Plotter([fits_file], save_path=output_dir)
        plot_fits.add_cphase().add_vis()
        # if mode == "incoherent":
            # plot_fits.add_flux()
        plot_fits.plot(save=True)

    cprint(f"{'':-^50}", "lg")
    cprint("Done!", "g")
    cprint(f"{'':-^50}", "lg")


def calibrate_folders(root_dir: Path, band_dir: Path,
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
    sub_dirs = sorted((root_dir / band_dir).glob("*.rb"))
    sub_dirs_rotated = deque(sub_dirs.copy())
    sub_dirs_rotated.rotate(1)

    for directory in sub_dirs:
        cprint(f"Calibration of {directory.name} with mode_name={mode_name}", "lp")
        cprint(f"{'':-^50}", "lg")
        if check_if_target(directory):
            for dir_rotated in sub_dirs_rotated:
                calibrate_fits_files(root_dir, directory,
                                     dir_rotated, mode_name=mode_name)
        else:
            cprint("No 'TARGET_RAW_INT*'-files found. SKIPPED!", "y")
            cprint(f"{'':-^50}", "lg")
            continue


def calibrate(data_dir: Path, stem_dir: Path, target_dir: Path):
    """Does the full calibration for all of the "cal_dir" subdirectories

    Parameters
    ----------
    data_dir: Path
    stem_dir: Path
    target_dir: Path
    """
    root_dir = Path(data_dir, stem_dir, "products", target_dir)
    modes, bands = {"coherent": "corrflux", "incoherent": "flux"}, ["lband", "nband"]

    for mode, mode_name in modes.items():
        for band in bands:
            calibrate_folders(root_dir, Path(mode, band), mode_name=mode_name)
    cprint("Calibration Done!", "lp")


if __name__ == "__main__":
    data_dir = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"
    stem_dir, target_dir = "matisse/GTO/hd163296/", "ATs/20190323"
    calibrate(data_dir, stem_dir, target_dir)
