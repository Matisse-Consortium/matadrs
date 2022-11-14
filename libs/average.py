import os
import matplotlib.pyplot as plt

from glob import glob
from pathlib import Path
from typing import Optional

from calib_BCD2 import calib_BCD
from avg_oifits import avg_oifits
from plot import Plotter

# NOTE: non chopped exposures 5-6. Do not average those together (For L-band)

def single_average(root_dir: Path, mode: str,
                   do_plots: Optional[bool] = False) -> None:
        """Calls Jozsef's code and does a average over the files for one band to average
        the reduced and calibrated data

        root_dir: Path
            The root folder for the PRODUCT
        band_dir: Path
            The folder in which the calibrated data is contained
        do_plots: bool, optional
        """
        # band_name = "L" if "lband" in band_dir else "N"
        mode_dir = os.path.join("calib", mode)

        folders = glob(os.path.join(root_dir, mode_dir , "*.rb"))
        for folder in folders:
            print(f"Averaging folder {os.path.basename(folder)}")
            fits_files = glob(os.path.join(folder, "*.fits"))
            fits_files.sort(key=lambda x: x[-8:])
            fits_files = fits_files[:4]

            folder_split = os.path.basename(folder).split(".")
            folder_split[0] += "-AVG"
            new_folder = ".".join(folder_split)
            outfile_dir = os.path.join(root_dir, "bcd_and_averaged",
                                       mode, new_folder)

            if not os.path.exists(outfile_dir):
                os.makedirs(outfile_dir)

            # epoch = ".".join(os.path.basename(band_dir).split(".")[2:5:2])
            outfile_path_vis = os.path.join(outfile_dir, "TARGET_AVG_VIS_INT.fits")
            outfile_path_cphases = os.path.join(outfile_dir, "TARGET_AVG_T3PHI_INT.fits")
            avg_oifits(fits_files, outfile_path_vis)
            calib_BCD(*fits_files, outfile_path_cphases, plot=False)

            # TODO: Make this with better plotter functionality
            if do_plot:
                print("Creating plots...")
                print("------------------------------------------------------------")
                fig = plt.figure()
                ax = fig.add_subplot()
                plot_vis = Plotter(outfile_path_vis).plot_vis(ax)
                plt.savfig("TARGET_AVG_VIS_INT.pdf", format="pdf")

                fig = plt.figure()
                ax = fig.add_subplot()
                plot_cphases = Plotter(outfile_path_cphases).plot_cphases(ax)
                plt.savfig("TARGET_AVG_T3PHI_INT.pdf", format="pdf")
            print("Done!")
            print("------------------------------------------------------------")


def averaging_pipeline(data_path: Path, stem_dir: Path,
                       target_dir: Path, do_plots: Optional[bool] = False) -> None:
        """Calls Jozsef's code sequentially to average the reduced and calibrated data

        Parameters
        ----------
        root_dir: Path
            The folder to the time-stamp folder
        do_plots: bool, optional
        """
        root_dir = os.path.join(data_path, stem_dir, "PRODUCT", target_dir)
        for mode in ["coherent", "incoherent"]:
            print(f"Averaging and BCD-calibration of {stem_dir} with mode={mode}")
            print("-----------------------------------------------------------")
            single_average(root_dir, mode)


if __name__ == "__main__":
    data_path = "/data/beegfs/astro-storage/groups/matisse/scheuck/data"
    stem_dir, target_dir = "matisse/GTO/hd142666/", "UTs/20220422"
    averaging_pipeline(data_path, stem_dir, target_dir)

