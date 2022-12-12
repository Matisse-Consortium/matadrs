from typing import List
from pathlib import Path
from collections import namedtuple

from .calib_BCD2 import calib_BCD
from .avg_oifits import avg_oifits

from .plot import Plotter
from .readout import ReadoutFits
from .utils import cprint, oifits_patchwork


# TODO: Make the files more modular
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


def merge_vis_and_cphases(stem_dir: Path, average_dir: Path) -> str:
    """Merges the vis and cphases files in the respective directory"""
    target_name = stem_dir.split("/")[~1]
    epoch = average_dir.name.split(".")[2]
    lband = True if "HAWAII" in average_dir else False
    band = "L" if lband else "N"
    fits_files = average_dir.glob("*.fits")
    cphases_file = [directory for directory in fits_files if "t3" in directory.lower()][0]
    vis_file  = [directory for directory in fits_files if "vis" in directory.lower()][0]
    out_file = average_dir / f"{target_name}_{epoch}_{band}_TARGET_AVG_INT.fits"
    oifits_patchwork(cphases_file, vis_file, out_file)
    plot = Plotter([out_file], lband=lband, save_path=out_file.parent)
    plot.add_cphases().add_corr_flux().plot(save=True)
    return out_file


def single_average(root_dir: Path, mode: str, lband) -> None:
        """Calls Jozsef's code and does a average over the files for one band to average
        the reduced and calibrated data

        root_dir: Path
            The root folder for the PRODUCT
        band_dir: Path
            The folder in which the calibrated data is contained
        """
        mode_dir = Path("calib", mode)
        folders = (root_dir / mode_dir).glob("*.rb")

        for folder in folders:
            print(f"Averaging folder {folder.name}")
            lband = True if "HAWAII" in folder else False

            unchopped_fits = folder.glob("*.fits")
            unchopped_fits.sort(key=lambda x: x[-8:])

            if len(fits_files) == 6:
                unchopped_fits = unchopped_fits[:4]
                chopped_fits = unchopped_fits[4:]

            folder_split = folder.name.split(".")
            folder_split[0] += "-AVG"
            new_folder = ".".join(folder_split)
            outfile_dir = root_dir / "bcd_and_averaged" / mode / new_folder

            if not outfile_dir.exists():
                outfile_dir.mkdir()

            outfile_path_vis = outfile_dir/ "TARGET_AVG_VIS_INT.fits"
            outfile_path_cphases = outfile_dir / "TARGET_AVG_T3PHI_INT.fits"
            outfile_path_chopped = outfile_dir / "TARGET_AVG_VIS_INT_CHOPPED.fits"
            bcd = sort_fits_by_BCD(unchopped_fits)
            avg_oifits(unchopped_fits, outfile_path_vis)

            if chopped_fits:
                avg_oifits(chopped_fits, outfile_path_chopped)

            if lband:
                calib_BCD(bcd.in_in, bcd.in_out, bcd.out_in, bcd.out_out,
                          outfile_path_cphases, plot=False)
            else:
                calib_BCD(bcd.in_in, "", "", bcd.out_out,
                         outfile_path_cphases, plot=False)

            # TODO: Make the plotter take the save name automatically,
            # if none given
            print("Creating plots...")
            print("------------------------------------------------------------")
            plot_vis = Plotter([outfile_path_vis], save_path=outfile_dir,
                               lband=lband, plot_name="TARGET_AVG_VIS_INT.pdf")
            plot_vis.add_vis().plot(save=True)

            plot_cphases = Plotter([outfile_path_cphases],
                                   save_path=outfile_dir, lband=lband,
                                   plot_name= "TARGET_AVG_T3PHI_INT.pdf")
            plot_cphases.add_cphases().plot(save=True)
            print("Plots created!")
            print("------------------------------------------------------------")
            print("Done!")
            print("-----------------------------------------------------------")


def average(data_path: Path, stem_dir: Path, target_dir: Path):
        """Calls Jozsef's code sequentially to average the reduced and calibrated data

        Parameters
        ----------
        root_dir: Path
            The folder to the time-stamp folder
        """
        root_dir = Path(data_path, stem_dir, "products", target_dir)
        for mode in ["coherent", "incoherent"]:
            cprint(f"Averaging and BCD-calibration of {stem_dir} with mode={mode}",
                   "lp")
            cprint("-----------------------------------------------------------",
                   "lg")
            single_average(root_dir, mode)
        cprint("Averaging done!", "lp")
        cprint("-----------------------------------------------------------",
               "lg")


if __name__ == "__main__":
    data_path = "/data/beegfs/astro-storage/groups/matisse/scheuck/data"
    stem_dir, target_dir = "matisse/GTO/hd163296/", "ATs/20190323"
    average(data_path, stem_dir, target_dir)

