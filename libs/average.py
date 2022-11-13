import os

from glob import glob
from pathlib import Path
from typing import List, Optional

from calib_BCD2 import calib_BCD
from avg_oifits import avg_oifits

# The data path to the general data
DATA_PATH = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"

# TODO: Comment and understand this part of Jozsef's code


def single_average(calib_dir: Path, bands: List[str]) -> None:
        """Calls Jozsef's code and does a average over the files for one band to average
        the reduced and calibrated data

        calib_dir: Path
            The folder in which the calibrated data is contained
        bands: List[str]
            A list containing the bands to be averaged over
        """
        for band in bands:
            band_dirs = glob(os.path.join(calib_dir, f"combined/{band}", "*.rb"))
            outfile_dir = os.path.join(os.path.basename(calib_dir),
                                       "combined_and_averaged")

            if not os.path.exists(outfile_dir):
                os.makedirs(outfile_dir)

            band_name = "L" if band == "lband" else "N"

            for band_dir in band_dirs:
                print(f"Averaging all files in {os.path.basename(band_dir)}")
                print("------------------------------------------------------------")
                fits_files = glob(os.path.join(band_dir, "*.fits"))
                fits_files.sort(key=lambda x: x[-8:])

                epoch = ".".join(os.path.basename(band_dir).split(".")[2:5:2])
                outfile_path_vis = os.path.join(outfile_dir,
                                                f"{epoch}_cal_avg_{band_name}_vis.fits")
                outfile_path_cphases = os.path.join(outfile_dir,
                                                    f"{epoch}_cal_avg_{band_name}_cphases.fits")
                avg_oifits(fits_files, outfile_path_vis)
                calib_BCD(*fits_files, outfile_path_cphases)
                print("Done!")
                print("------------------------------------------------------------")

def averaging_pipeline(calib_dir: Path,
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
        # non chopped exposures 5-6. Do not average those together
        if both:
            bands = ["lband", "nband"]
        else:
            bands = ["lband"] if lband else ["nband"]

        single_average(calib_dir, bands)


if __name__ == "__main__":
    specific_dir = "GTO/hd142666/PRODUCT/UTs/20190514"
    averaging_pipeline(os.path.join(DATA_PATH, specific_dir), both=True)

