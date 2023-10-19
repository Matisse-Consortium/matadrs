import re
import shutil
import subprocess
from collections import deque, namedtuple
from pathlib import Path
from typing import Optional, Dict, List
from warnings import warn

import pkg_resources

from .fluxcal import fluxcal
from .calib_BCD2 import calib_BCD
from ..utils.plot import Plotter
from ..utils.readout import ReadoutFits
from ..utils.tools import cprint, print_execution_time, get_path_descriptor,\
        check_if_target, get_fits_by_tag, get_execution_modes, split_fits

__all__ = ["create_visibility_sof", "check_file_match", "sort_fits_by_bcd_configuration",
           "calibrate_bcd", "calibrate_visibilities", "calibrate_fluxes",
           "cleanup_calibration", "calibrate_files", "calibrate_folders", "calibrate"]

DATABASE_DIR = Path(pkg_resources.resource_filename("matadrs", "data/calibrator_databases"))
DATABASES = ["vBoekelDatabase.fits", "calib_spec_db_v10.fits",
             "calib_spec_db_v10_supplement.fits", "calib_spec_db_supplement3.fits"]
LBAND_DATABASES = list(map(lambda x: DATABASE_DIR / x, DATABASES))
NBAND_DATABASES = LBAND_DATABASES[:]+[DATABASE_DIR / "vBoekelDatabase.fitsold"]

MODE_NAMES = {"coherent": "corrflux", "incoherent": "flux"}


# TODO: Get correct calibrator for N- and L-band
# TODO: Make functionality that does not only calibrate willy nilly, but checks if the
# calibrator is there for N-band, L-band or LN-band, see Jozsef's files. Use the
# mat_target_list of newest edition for that (grep it via python api of google sheets)
# Make a function for this
# MAT_TARGET_LIST = DATA_DIR / "mat_target_list.xlsx"
def create_visibility_sof(reduced_dir: Path,
                          targets: List[Path],
                          calibrators: List[Path]) -> Path:
    """Creates the (.sof)-file needed for the visibility calibration with
    'mat_cal_oifits' and returns its path.

    Parameters
    ----------
    reduced_dir : pathlib.Path
        The directory containing the, by the pipeline, reduced data.
    targets : list of pathlib.Path
        The 'TARGET_RAW_INT'-files.
    calibrators : list of pathlib.Path
        The 'CALIB_RAW_INT'-files.

    Returns
    -------
    sof_file : pathlib.Path
        The Path to the (.sof)-file used for visibility calibration.
    """
    sof_file = reduced_dir / "visbility_reduction.sof"
    with open(sof_file, "w+") as sof:
        for target in targets:
            sof.write(f"{target} TARGET_RAW_INT\n")
        sof.write("\n")
        for calibrator in calibrators:
            sof.write(f"{calibrator} CALIB_RAW_INT\n")
    return sof_file


def check_file_match(targets: List[Path], calibrators: List[Path]) -> bool:
    """Checks if all needed files exist in the directories and if the number
    of both directories is the same.

    Parameters
    ----------
    targets : list of pathlib.Path
        The detected "TARGET_RAW"-files.
    calibrators : list of pathlib.Path
        The detected "CALIB_RAW"-files.

    Returns
    -------
    files_match : bool
    """
    if not targets:
        cprint("No 'TARGET_RAW_INT*'-files found (Maybe calibrator? If not check for"
               " error in first reduction step). SKIPPING!", "y")
        cprint(f"{'':-^50}", "lg")
        return False
    if len(targets) < 4:
        warn("# 'TARGET_RAW_INT'-files is lower than 4! Indicates problems with reduction. SKIPPING!")
        return False
    if len(targets) != len(calibrators):
        warn("#'TARGET_RAW_INT'-files != #'CALIB_RAW_INT'-files!")
    return True


def sort_fits_by_bcd_configuration(fits_files: List[Path]) -> namedtuple:
    """Sorts the input (.fits)-files by their BCD configuration.

    Parameters
    ----------
    fits_files : list of pathlib.Path
        The (.fits)-files to be checked for their BCD-configuration.

    Returns
    -------
    bcd_configuration : collections.namedtuple
        A named tuple containing the (.fits)-files sorted by BCD-configuration.
        If none of the specified configuration is found then the field is
        filled with an empty string.
    """
    in_in, in_out, out_in, out_out = "", "", "", ""
    BCDFits = namedtuple("BCDFits", ["in_in", "in_out", "out_in", "out_out"])
    for fits_file in fits_files:
        bcd_configuration = ReadoutFits(fits_file).bcd_configuration
        if bcd_configuration == "in-in":
            in_in = fits_file
        elif bcd_configuration == "in-out":
            in_out = fits_file
        elif bcd_configuration == "out-in":
            out_in = fits_file
        elif bcd_configuration == "out-out":
            out_out = fits_file
        else:
            cprint("BCD-configuration has not been found!", "r")
            raise ValueError
    return BCDFits(in_in, in_out, out_in, out_out)


def match_targets_to_calibrators(targets: List[Path],
                                 calibrators: List[Path]) -> Dict[str, str]:
    """Matches the 'TARGET_RAW_INT'- to the 'CALIB_RAW_INT'-files."""
    first_dict = {}
    for target in targets:
        numerical_part = re.search(r'\d+', target.name).group()
        first_dict[numerical_part] = target

    matched_entries = {}
    for calibrator in calibrators:
        numerical_part = re.search(r'\d+', calibrator.name).group()
        if numerical_part in first_dict:
            matched_entries[first_dict[numerical_part]] = calibrator
    return matched_entries


def calibrate_bcd(directory: Path, band: str, output_dir: Path) -> None:
    """Executes the BCD-calibration for the the unchopped/chopped, visbility calibrated
    files.

    Parameters
    ----------
    directory : pathlib.Path
        The directory to be searched in.
    output_dir : pathlib.Path
        The directory to which the new files are saved to.

    Notes
    -----
    This creates either one or two output files depending if there are only
    unchopped or also chopped files. The files' names end with either 'INT' or
    'INT_CHOPPED', respectively, and indicated that they are averaged by an
    'AVG' in their name.

    From MATISSE-pipeline version 1.7.6, the BCD-calibration is included in the
    pipeline, however, for L-band it only takes into account the chopped files,
    and is unreliable if they are not present, in which case the
    BCD-calibration is executed via this script and the output of it merged
    into the final file.

    See also
    --------
    .calib_BCD2.calib_BCD : BCD-calibration for closure phases.
    """
    cprint("Executing BCD-calibration...", "g")
    unchopped_fits, chopped_fits = split_fits(directory, "CAL_INT_0")
    outfile_unchopped_cphases = output_dir / "TARGET_BCD_CAL_T3PHI_INT.fits"
    bcd = sort_fits_by_bcd_configuration(unchopped_fits)
    if band == "lband":
        calib_BCD(bcd.in_in, bcd.in_out,
                  bcd.out_in, bcd.out_out,
                  outfile_unchopped_cphases, plot=False)

        if chopped_fits is not None:
            outfile_chopped_cphases = output_dir / "TARGET_BCD_CAL_T3PHI_INT_CHOPPED.fits"
            bcd_chopped = sort_fits_by_bcd_configuration(chopped_fits)
            calib_BCD(bcd_chopped.in_in, "",
                      "", bcd_chopped.out_out,
                      outfile_chopped_cphases, plot=False)
    else:
        calib_BCD(bcd.in_in, "", "", bcd.out_out, outfile_unchopped_cphases, plot=False)


def calibrate_visibilities(targets: List[Path],
                           calibrators: List[Path], output_dir: Path) -> None:
    """Calibrates the visibilities of all the provided files and saves them to
    the output directory.

    Parameters
    ----------
    targets : list of pathlib.Path
        The detected "TARGET_RAW"-files.
    calibrators : list of pathlib.Path
        The detected "CALIB_RAW"-files.
    product_dir : pathlib.Path
        The directory to contain the calibrated files.
    """
    cprint("Calibrating visibilities...", "g")
    sof_file = create_visibility_sof(output_dir, targets, calibrators)
    subprocess.call(["esorex", f"--output-dir={str(output_dir)}",
                     "mat_cal_oifits", str(sof_file)],
                    stdout=subprocess.DEVNULL)
    cprint("Plotting visibility calibrated files...", "y")
    for fits_file in get_fits_by_tag(output_dir, "TARGET_CAL_INT"):
        plot_fits = Plotter(fits_file, save_path=output_dir)
        plot_fits.add_cphases().add_vis().plot(save=True, error=True)


def calibrate_fluxes(targets: List[Path], calibrators: List[Path],
                     mode: str, band: str, output_dir: Path) -> None:
    """Calibrates the fluxes of all the provided files and saves it to the
    output directory.

    Parameters
    ----------
    targets : list of pathlib.Path
        The detected "TARGET_RAW"-files.
    calibrators : list of pathlib.Path
        The detected "CALIB_RAW"-files.
    mode : str
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'.
    band : str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'.
    product_dir : pathlib.Path
        The directory to contain the calibrated files.
    """
    cprint("Calibrating fluxes...", "g")
    for index, (target, calibrator) in enumerate(zip(targets, calibrators), start=1):
        cprint(f"Processing {target.name} with {calibrator.name}...", "g")
        output_file = output_dir / f"TARGET_FLUX_CAL_INT_000{index}.fits"
        databases, do_airmass = NBAND_DATABASES, True
        if band == "lband":
            databases, do_airmass = LBAND_DATABASES, False
        fluxcal(str(target), str(calibrator), str(output_file),
                list(map(str, databases)), mode=MODE_NAMES[mode],
                output_fig_dir=str(output_dir),
                do_airmass_correction=do_airmass)
        cprint(f"Plotting file '{output_file.name}'...", "y")
        plot_fits = Plotter(output_file, save_path=output_dir)
        plot_fits.add_cphases().add_vis(corr_flux=True).add_vis2().plot(save=True, error=True)


# TODO: Find way to make this moving better than this -> Moves (.fits)-files
# from this path that get created by 'mat_cal_oifits?'
def cleanup_calibration(output_dir: Path):
    """Moves any (.fits)-files created by 'esorex' and the esorex (.log)-files
    to the calibrated folder."""
    for fits_file in Path().cwd().glob("*.fits"):
        shutil.move(str(fits_file), str(output_dir / fits_file.name))
    shutil.move(str(Path().cwd() / "esorex.log"),
                str(output_dir / "mat_cal_oifits.log"))
    cprint(f"{'':-^50}", "lg")


def calibrate_files(reduced_dir: Path, target_dir: Path,
                    calibrator_dir: Path, mode: str,
                    band: str, overwrite: bool) -> None:
    """The total calibration for a target and a calibrator folder. Includes
    the flux-, visibility- and closure phase (bcd-) calibration.

    Parameters
    ----------
    reduced_dir : pathlib.Path
        The path to multiple folders that need to be cross correlated.
    target_dir : pathlib.Path
        The directory of the science target.
    calibrator_dir : pathlib.Path
        The directory of the calibrator.
    mode : str
        The mode of calibration. Either 'corrflux', 'flux' or 'both' depending
        if it is 'coherent' or 'incoherent' reduced data.
    band : str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'.
    overwrite : bool, optional
        If 'True' overwrites files from previous calibration.
    """
    cprint(f"Calibrating {target_dir.name} with {calibrator_dir.name}...", "p")
    targets = get_fits_by_tag(target_dir, "TARGET_RAW_INT")
    calibrators = get_fits_by_tag(calibrator_dir, "CALIB_RAW_INT")

    if check_file_match(targets, calibrators):
        matches = match_targets_to_calibrators(targets, calibrators)
        targets, calibrators = map(list, [matches.keys(), matches.values()])
        output_dir = get_path_descriptor(reduced_dir, "TAR-CAL",
                                         targets[0], calibrators[0])
        if not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=overwrite)
        calibrate_fluxes(targets, calibrators, mode, band, output_dir)
        calibrate_visibilities(targets, calibrators, output_dir)
        calibrate_bcd(output_dir, band, output_dir)
        cleanup_calibration(output_dir)


def calibrate_folders(reduced_dir: Path, mode: str,
                      band: str, overwrite: bool) -> None:
    """Calibrates a directory containing the scientific target with a directory
    containing the calibrator observation. Calibrates flux, visibility and
    closure phases (bcd).

    Parameters
    ----------
    reduced_dir : pathlib.Path
        The path to multiple folders that need to be cross correlated.
    mode : str, optional
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'.
    band : str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'.
    overwrite : bool, optional
        If 'True' overwrites files from previous calibration.
    """
    sub_dirs = sorted((reduced_dir / "reduced" / mode / band).glob("*.rb"))
    rotated_sub_directories = deque(sub_dirs.copy())
    rotated_sub_directories.rotate(1)

    for directory in sub_dirs:
        cprint(f"Calibration of {directory.name} in '{mode}' mode", "lp")
        cprint(f"{'':-^50}", "lg")
        if check_if_target(directory):
            for rotated_directory in rotated_sub_directories:
                if directory == rotated_directory:
                    continue
                calibrate_files(reduced_dir, directory,
                                rotated_directory, mode, band, overwrite)
    cprint(f"Finished calibration of {band} and {mode}", "lp")
    cprint(f"{'':-^50}", "lp")


# TODO: Remove redundant checks for science targets as well as SCI-SCI or CAL-CAL. Skip
# cals as checker generally
# TODO: Implement checking for overwriting. Right now overwriting is by default
@print_execution_time
def calibration_pipeline(reduced_dir: Path,
                         mode: Optional[str] = "both",
                         band: Optional[str] = "both",
                         overwrite: Optional[bool] = False) -> None:
    """Does the full calibration for all of the reduced directories
    subdirectories.

    Parameters
    ----------
    reduced_dir : pathlib.Path
        The path to multiple folders that need to be cross correlated.
    mode : str, optional
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'.
    band : str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'.
    overwrite : bool, optional
        If 'True' overwrites files from previous calibration.

    Notes
    -----
    The mode of calibration can either be 'coherent' or 'incoherent', which
    corresponds to either the calibration of the correlated fluxes or the
    total flux.
    """
    modes, bands = get_execution_modes(mode, band)
    for mode in modes:
        for band in bands:
            calibrate_folders(reduced_dir, mode, band, overwrite)
    cprint(f"Finished calibration of {', '.join(bands)} and {', '.join(modes)}", "lp")


if __name__ == "__main__":
    data_dir = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"
    stem_dir, target_dir = "matisse/GTO/hd142666/", "UTs/20220420"
    calibrate(data_dir, stem_dir, target_dir)
