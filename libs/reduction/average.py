import shutil
from typing import List
from pathlib import Path
from collections import namedtuple

# TODO: Find way to make this into a complete module -> More pythonic!
from readout import ReadoutFits
from calib_BCD2 import calib_BCD
from avg_oifits import avg_oifits
from utils import cprint, get_fits


HEADER_TO_REMOVE = [{'key':'HIERARCH ESO INS BCD1 ID','value':' '},
                    {'key':'HIERARCH ESO INS BCD2 ID','value':' '},
                    {'key':'HIERARCH ESO INS BCD1 NAME','value':' '},
                    {'key':'HIERARCH ESO INS BCD2 NAME','value':' '}]


def sort_fits_by_BCD(fits_files: List[Path]) -> namedtuple:
    """Sorts the input (.fits)-files by their BCD configuration

    Parameters
    ----------
    fits_files: List[Path]

    Returns
    -------
    namedtuple
        The fits ordered by their BCD-configuration
    """
    BCDFits = namedtuple("BCDFits", ["in_in", "in_out", "out_in", "out_out"])
    # TODO: Can be easier done with DataPrep -> Implement here
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


def average_files(unchopped_fits: List[Path], chopped_fits: List[Path], output_dir: Path):
    """Averages a set of unchopped and chopped (.fits)-files"""
    if "TARGET_CAL_INT" in unchopped_fits[0].name:
        cprint(f"Averaging visibility calibration...", "g")
        outfile_name = "TARGET_AVG_VIS"
    else:
        cprint(f"Averaging flux calibration...", "g")
        outfile_name = "TARGET_AVG_FLUX"

    outfile_unchopped = output_dir/ f"{outfile_name}_INT.fits"
    outfile_chopped = output_dir/ f"{outfile_name}_INT_CHOPPED.fits"

    avg_oifits(unchopped_fits, outfile_unchopped, headerval=HEADER_TO_REMOVE)

    if chopped_fits is not None:
        avg_oifits(chopped_fits, outfile_chopped, headerval=HEADER_TO_REMOVE)


def bcd_calibration(unchopped_fits: List[Path], output_dir: Path):
    """Executes the BCD-calibration for the the unchopped visbility calibrated files"""
    cprint(f"Executing BCD-calibration...", "g")
    outfile_path_cphases = output_dir / "TARGET_AVG_T3PHI_INT.fits"
    bcd = sort_fits_by_BCD(unchopped_fits)
    if "lband" in str(unchopped_fits[0]):
        calib_BCD(bcd.in_in, bcd.in_out, bcd.out_in, bcd.out_out,
                  outfile_path_cphases, plot=False)
    else:
        calib_BCD(bcd.in_in, "", "", bcd.out_out,
                 outfile_path_cphases, plot=False)
    # TODO: See how to bcd-calibrate the chopped files as well


def average_folders(root_dir: Path, mode: str) -> None:
        """Calls Jozsef's code and does a average over the files for one band to average
        the reduced and calibrated data

        Parameters
        ----------
        root_dir: Path
            The root folder for the PRODUCT
        band_dir: Path
            The folder in which the calibrated data is contained
        """
        mode_dir = Path("calib", mode)
        folders = (root_dir / mode_dir).glob("*.rb")

        for folder in folders:
            # TODO: Make this into function
            cprint(f"Averaging folder {folder.name}...", "g")
            cprint(f"{'':-^50}", "lg")

            # TODO: Make check if vis-calibration or flux calibration did not work
            unchopped_vis_fits, chopped_vis_fits = get_fits(folder, "TARGET_CAL")
            unchopped_flux_fits, chopped_flux_fits = get_fits(folder, "TARGET_FLUXCAL")

            folder_split = folder.name.split(".")
            folder_split[0] += "-AVG"
            new_folder = ".".join(folder_split)
            output_dir = root_dir / "bcd_and_averaged" / mode / new_folder

            if not output_dir.exists():
                output_dir.mkdir(parents=True)

            average_files(unchopped_vis_fits, chopped_vis_fits, output_dir)
            average_files(unchopped_flux_fits, chopped_flux_fits, output_dir)
            bcd_calibration(unchopped_vis_fits, output_dir)
            for fits_file in folder.glob("*.fits"):
                shutil.copy(str(fits_file), (output_dir / fits_file.name))
            cprint(f"{'':-^50}", "lg")


def average(data_path: Path, stem_dir: Path, target_dir: Path):
        """Calls Jozsef's code sequentially to average the reduced and calibrated data

        Parameters
        ----------
        root_dir: Path
            The folder to the time-stamp folder
        """
        root_dir = Path(data_path, stem_dir, "products", target_dir)
        for mode in ["coherent", "incoherent"]:
            cprint(f"Averaging and BCD-calibration of {stem_dir} with mode={mode}", "lp")
            cprint(f"{'':-^50}", "lg")
            average_folders(root_dir, mode)
        cprint("Averaging done!", "lp")
        cprint(f"{'':-^50}", "lg")


if __name__ == "__main__":
    data_path = "/data/beegfs/astro-storage/groups/matisse/scheuck/data"
    stem_dir, target_dir = "matisse/GTO/hd163296/", "ATs/20190323"
    average(data_path, stem_dir, target_dir)
