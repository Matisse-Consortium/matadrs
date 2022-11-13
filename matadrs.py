import os
import time
import datetime

from pathlib import Path
from typing import List, Optional

from libs.reduction import reduction_pipeline
from libs.calibration import calibration_pipeline
from libs.average import averaging_pipeline
from libs.merge import merging_pipeline

# TODO: Add functionality that clears all the paths before it write again, as to overwrite them

# The data path to the general data
DATA_PATH = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"

def full_pipeline(raw_dirs: List[Path],
                  arrays: Optional[List[str]] = [],
                  res_dirs: Optional[List[Path]] = [],
                  calib_dirs: Optional[List[Path]] = [],
                  both_lst: Optional[List[bool]] = [],
                  lbands: Optional[List[bool]] = []) -> None:
    """Combines all the facettes of data reduction into one executable function that takes
    a list of epochs and different datasets to be reduced via the MATISSE pipeline, then
    calibrated, merged and averaged

    Parameters
    ----------
    raw_dirs: List[Path]
        The respective RAW directories
    arrays: List[str], optional
        A list of arrays that are used for the reduction. If left empty, then this function
        tries to find the correct configuration via the path that is given in the "raw_dirs".
        If nothing is found an error is raised
    res_dirs: List[Path], optional
        The respective PRODUCT directories. Will be set automatically if not given to
        "PRODUCTS" instead of raw, maintaining the same path as the "raw_dirs"
    calib_dirs: List[Path], optional
        The respective CALIB directories, will be set to the corresponding RAW directories
        if left empty
    both_lst: List[bool], optional
        A list containing the respective values for the reduction of the bands. If both
        bands are to be averaged, this has to be false for the "lband" option to be
        considered. If left empty then it is set to "True"
    lbands: List[bool], optional
        A list containting the respective values for the reduction of the "lband" or
        "nband". If "both=False" and this is "True"", then lband will be averaged over, if
        "both=False" and this is "False", then nband will be averaged over. If left empty
        then it is set to "False"

    See Also
    --------
    reduction_pipeline()
    calibration_pipeline()
    merging_pipeline()
    averaging_pipeline()
    """
    for i, raw_dir in enumerate(raw_dirs):
        raw_dir = os.path.join(DATA_PATH, raw_dir)

        if not arrays:
            if "ATs" in raw_dir:
                array = "ATs"
            elif "UTs" in raw_dir:
                array = "UTs"
            else:
                raise IOError("No array configuration has been input and it could"\
                              " also not be determined via the folder structure!")

        if not calib_dirs:
            calib_dir = raw_dir
        else:
            calib_dir = calib_dirs[i]

        if not res_dirs:
            res_dir = raw_dir.replace("RAW", "PRODUCT")
        else:
            res_dir = res_dirs[i]

        if not both_lst:
            both, lband = True, False
        else:
            both, lband = both_lst[i], lbands[i]

        start_time = time.time()
        print(f"Starting data improvement of {raw_dir}!")
        print("----------------------------------------------------------------------")
#        reduction_pipeline(raw_dir, calib_dir, res_dir, array,
#                           lband=lband, both=both)
        calibration_pipeline(res_dir, lband=lband, both=both)
        merging_pipeline(res_dir, lband=lband, both=both)
        averaging_pipeline(res_dir, lband=lband, both=both)
        print("----------------------------------------------------------------------")
        print(f"Reduction, calibration, merging and averaging complete in "\
              f"{datetime.timedelta(seconds=(time.time()-start_time))} hh:mm:ss")
    print("----------------------------------------------------------------------")
    print("All done!")


if __name__ == "__main__":
    specific_path = "GTO/hd142666/RAW/"
    uts = [os.path.join("UTs", i) for i in ["20190514", "20220420", "20220422"]]
#    ats = ["medium/20190506", "medium/20220531",
#           "medium/20220601", "small/20190323", "small/20190630"]
#    ats = [os.path.join("ATs", i) for i in ats]
    raw_dirs = []
    raw_dirs.extend(uts)
    raw_dirs = [os.path.join(specific_path, i) for i in raw_dirs]
    full_pipeline(raw_dirs)
