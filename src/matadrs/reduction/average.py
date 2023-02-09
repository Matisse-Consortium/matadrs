"""

Routine
-------

See Also
--------

References
----------

Examples
--------
"""

__all__ = ["sort_fits_by_bcd_configuration", "average_files",
           "bcd_calibration", "average_folders", "average"]

import shutil
from collections import namedtuple
from pathlib import Path
from typing import Optional, List

from .avg_oifits import avg_oifits
from .calib_BCD2 import calib_BCD
from ..utils.plot import Plotter
from ..utils.readout import ReadoutFits
from ..utils.tools import cprint, split_fits, get_fits_by_tag, get_execution_modes


HEADER_TO_REMOVE = [{'key': 'HIERARCH ESO INS BCD1 ID', 'value': ' '},
                    {'key': 'HIERARCH ESO INS BCD2 ID', 'value': ' '},
                    {'key': 'HIERARCH ESO INS BCD1 NAME', 'value': ' '},
                    {'key': 'HIERARCH ESO INS BCD2 NAME', 'value': ' '}]


def sort_fits_by_bcd_configuration(fits_files: List[Path]) -> namedtuple:
    """Sorts the input (.fits)-files by their BCD configuration

    Parameters
    ----------
    fits_files: List[Path]
        The (.fits)-files to be checked for their BCD-configuration

    Returns
    -------
    bcd_configuration: namedtuple
        A named tuple containing the (.fits)-files sorted by BCD-configuration. If none of
        the specified configuration is found then the field is filled with an empty string
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


def average_files(directory: Path, file_type: str, output_dir: Path) -> None:
    """Averages the unchopped (.fits)-files and the chopped (.fits)-files if given

    Parameters
    ----------
    directory: Path
        The directory to be searched in
    file_type: str
        Either 'flux' or 'vis' for the flux and visibility calibration, respectively
    output_dir: Path
        The directory to which the new files are saved to

    Notes
    -----
    This creates either one or two output files depending if there are only unchopped or
    also chopped files. The files' names are with either 'INT' or 'INT_CHOPPED',
    respectively, and indicated that they are averaged by an 'AVG' in their name


    See also
    --------
    .avg_oifits.avg_oifits: Averaging for (.fits)-files
    """
    if file_type == "flux":
        cprint("Averaging flux calibration...", "g")
        outfile_name = "TARGET_AVG_FLUX"
        unchopped_fits, chopped_fits = split_fits(directory, "TARGET_FLUX_CAL")
    else:
        cprint("Averaging visibility calibration...", "g")
        unchopped_fits, chopped_fits = split_fits(directory, "TARGET_CAL")
        outfile_name = "TARGET_AVG_VIS"

    outfile_unchopped = output_dir / f"{outfile_name}_INT.fits"
    avg_oifits(unchopped_fits, outfile_unchopped, headerval=HEADER_TO_REMOVE)

    if chopped_fits is not None:
        outfile_chopped = output_dir / f"{outfile_name}_INT_CHOPPED.fits"
        avg_oifits(chopped_fits, outfile_chopped, headerval=HEADER_TO_REMOVE)


def bcd_calibration(directory: Path, output_dir: Path) -> None:
    """Executes the BCD-calibration for the the unchopped, visbility calibrated files

    Parameters
    ----------
    directory: Path
        The directory to be searched in
    output_dir: Path
        The directory to which the new files are saved to

    Notes
    -----
    This creates either one or two output files depending if there are only unchopped or
    also chopped files. The files' names end with either 'INT' or 'INT_CHOPPED',
    respectively, and indicated that they are averaged by an 'AVG' in their name

    See also
    --------
    .calib_BCD2.calib_BCD: BCD-calibration for closure phases
    """
    cprint("Executing BCD-calibration...", "g")
    unchopped_fits, chopped_fits = split_fits(directory, "TARGET_CAL")
    outfile_unchopped_cphases = output_dir / "TARGET_AVG_T3PHI_INT.fits"
    bcd = sort_fits_by_bcd_configuration(unchopped_fits)
    if "lband" in str(unchopped_fits[0]):
        calib_BCD(bcd.in_in, bcd.in_out,
                  bcd.out_in, bcd.out_out,
                  outfile_unchopped_cphases, plot=False)
    else:
        calib_BCD(bcd.in_in, "", "", bcd.out_out, outfile_unchopped_cphases, plot=False)

    if chopped_fits is not None:
        outfile_chopped_cphases = output_dir / "TARGET_AVG_T3PHI_CHOPPED_INT.fits"
        bcd_chopped = sort_fits_by_bcd_configuration(chopped_fits)
        calib_BCD(bcd_chopped.in_in, "", "", bcd_chopped.out_out,
                  outfile_chopped_cphases, plot=False)


def average_folders(calibrated_dir: Path, mode: str) -> None:
        """Iterates over the calibrated directories to

        Parameters
        ----------
        calibrated_dir: Path
            The directory containing the calibration's products
        mode: str, optional
            The mode in which the reduction is to be executed. Either 'coherent',
            'incoherent' or 'both'
        """
        for directory in (calibrated_dir / "calib" / mode).glob("*.rb"):
            cprint(f"Averaging folder {directory.name}...", "g")
            cprint(f"{'':-^50}", "lg")

            folder_split = directory.name.split(".")
            folder_split[0] += "-AVG"
            new_folder = ".".join(folder_split)
            output_dir = calibrated_dir / "bcd_and_averaged" / mode / new_folder

            if not output_dir.exists():
                output_dir.mkdir(parents=True)

            average_files(directory, "flux", output_dir)
            average_files(directory, "vis", output_dir)
            bcd_calibration(directory, output_dir)
            for fits_file in get_fits_by_tag(directory, "CAL_INT"):
                shutil.copy(str(fits_file), (output_dir / fits_file.name))

            cprint("Plotting averaged files...", "g")
            for fits_file in get_fits_by_tag(output_dir, "AVG"):
                plot_fits = Plotter([fits_file], save_path=output_dir)
                plot_fits.add_cphases().add_vis().plot(save=True)
            cprint(f"{'':-^50}", "lg")


# TODO: Implement overwrite
def average(calibrated_dir: Path,
            mode: Optional[str] = "both", overwrite: Optional[bool] = False):
    """Does the full averaging for all of the calibrated directories subdirectories

    Parameters
    ----------
    calibrated_dir: Path
        The directory containing the calibration's products
    mode: str, optional
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'
    overwrite: bool, optional
        If 'True' overwrites present files from previous calibration
    """
    for mode in get_execution_modes(mode)[0]:
        cprint("Averaging and BCD-calibration of"
               f" {calibrated_dir.name} with mode={mode}", "lp")
        cprint(f"{'':-^50}", "lg")
        average_folders(calibrated_dir, mode)
    cprint("Averaging done!", "lp")
    cprint(f"{'':-^50}", "lg")
