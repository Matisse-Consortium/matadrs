from typing import List
from pathlib import Path
from collections import namedtuple

# TODO: Find way to make this into a complete module -> More pythonic!
from plot import Plotter
from readout import ReadoutFits
from calib_BCD2 import calib_BCD
from avg_oifits import avg_oifits
from utils import cprint, oifits_patchwork


HEADER_TO_REMOVE = [{'key':'HIERARCH ESO INS BCD1 ID','value':' '},
                    {'key':'HIERARCH ESO INS BCD2 ID','value':' '},
                    {'key':'HIERARCH ESO INS BCD1 NAME','value':' '},
                    {'key':'HIERARCH ESO INS BCD2 NAME','value':' '}]

def get_fits(folder: Path, tag: str):
    """Searches a folder for a tag and returns the non-chopped and chopped (.fits)-files"""
    unchopped_fits = folder.glob(f"{tag}*.fits")
    unchopped_fits.sort(key=lambda x: x[-8:])

    if len(unchopped_fits) == 6:
        return unchopped_fits[:4], unchopped_fits[4:]
    return unchopped_fits, None


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
    for fits_file in fits_files:
        readout = ReadoutFits(fits_file)
        bcd_info = readout.get_bcd_info()
        if bcd_info == "in-in":
            in_in = fits_file
        elif bcd_info == "in-out":
            in_out = fits_file
        elif bcd_info == "out-in":
            out_in = fits_file
        elif bcd_info == "out-out":
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

    # TODO: See how to bcd-calibrate the chopped files as well
    if chopped_fits is not None:
        avg_oifits(chopped_fits, outfile_chopped, headerval=HEADER_TO_REMOVE)
        # calib_BCD()


def bcd_calibration(unchopped_fits: List[Path], output_dir: Path):
    """Executes the BCD-calibration for the the unchopped visbility calibrated files"""
    cprint(f"Executing BCD-calibration...", "g")
    lband = True if "HAWAII" in unchopped_fits[0].parent else False
    outfile_path_cphases = output_dir / "TARGET_AVG_T3PHI_INT.fits"
    bcd = sort_fits_by_BCD(unchopped_fits)
    if lband:
        calib_BCD(bcd.in_in, bcd.in_out, bcd.out_in, bcd.out_out,
                  outfile_path_cphases, plot=False)
    else:
        calib_BCD(bcd.in_in, "", "", bcd.out_out,
                 outfile_path_cphases, plot=False)


def average_folders(root_dir: Path, mode: str, lband) -> None:
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

            # TODO: Make this work again
            # cprint("Creating plots...", "g")
            # cprint(f"{'':-^50}", "lg")
            # plot_vis = Plotter([outfile_path_vis], save_path=output_dir,)
            cprint(f"{'':-^50}", "lg")
            cprint("Done!", "g")
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
