import os

from glob import glob
from pathlib import Path
from typing import List, Optional

from calib_BCD2 import calib_BCD
from avg_oifits import avg_oifits

# The data path to the general data
DATA_PATH = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"


def single_average(root_dir: Path, band_dir: Path) -> None:
        """Calls Jozsef's code and does a average over the files for one band to average
        the reduced and calibrated data

        root_dir: Path
            The root folder for the PRODUCT
        band_dir: Path
            The folder in which the calibrated data is contained
        """
        band_name = "L" if "lband" in band_dir else "N"

        folders = glob(os.path.join(root_dir, band_dir , "*.rb"))
        for folder in folders:
            print(f"Averaging folder {os.path.basename(folder)}")
            fits_files = glob(os.path.join(folder, "*.fits"))
            fits_files.sort(key=lambda x: x[-8:])
            fits_files = fits_files[:4]

            folder_split = os.path.basename(folder).split(".")
            folder_split[0] += "-AVG"
            new_folder = ".".join(folder_split)
            outfile_dir = os.path.join(root_dir, "bcd_and_averaged", new_folder)

            if not os.path.exists(outfile_dir):
                os.makedirs(outfile_dir)

            epoch = ".".join(os.path.basename(band_dir).split(".")[2:5:2])
            outfile_path_vis = os.path.join(outfile_dir, "TARGET_AVG_VIS_INT.fits")
            outfile_path_cphases = os.path.join(outfile_dir, "TARGET_AVG_T3PHI_INT.fits")
            avg_oifits(fits_files, outfile_path_vis)
            calib_BCD(*fits_files, outfile_path_cphases, plot=False)
            print("Done!")
            print("------------------------------------------------------------")

def averaging_pipeline(root_dir: Path,
                       both: Optional[bool] = False,
                       lband: Optional[bool] = False) -> None:
        """Calls Jozsef's code sequentially to average the reduced and calibrated data

        Parameters
        ----------
        calib_dir: Path
            The folder in which the calibrated data is contained
        both: bool, optional
            If both bands are to be averaged, this has to be false for the "lband" option
            to be considered
        lband: bool, optional
            If "both=False" and this is "True"", then lband will be averaged over, if
            "both=False" and this is "False", then nband will be averaged over
        """
        modes, bands = ["coherent", "incoherent"], ["lband", "nband"]

        # non chopped exposures 5-6. Do not average those together (For L-band)
        if both:
            for mode in modes:
                mode_dir = os.path.join("calib", mode)
                for band in bands:
                    band_dir = os.path.join(mode_dir, band)
                    print(f"Averaging and BCD-calibration of {band_dir}")
                    print(f"with mode={mode}")
                    print("------------------------------------------------------------")
                    single_average(root_dir, band_dir)
        else:
            band = "lband" if lband else "nband"
            for mode in modes:
                band_dir = os.path.join("calib", mode, band)
                print(f"Averaging and BCD-calibration of {band_dir}")
                print(f"with mode={mode}")
                print("------------------------------------------------------------")
                single_average(root_dir, band_dir)


if __name__ == "__main__":
    specific_dir = "matisse/GTO/hd142666/PRODUCT/UTs/20220422"
    averaging_pipeline(os.path.join(DATA_PATH, specific_dir), both=True)

