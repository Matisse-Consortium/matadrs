import os
import time
import datetime

from pathlib import Path
from typing import Optional

from .reduction import reduction_pipeline, calibration_pipeline, averaging_pipeline,\
        merging_pipeline

# TODO: Add functionality that clears all the paths before it write again, as to overwrite them

def matadrs_pipeline(data_dir: Path, stem_dir: Path,
                     target_dirs: List[Path],
                     do_reduce: Optional[bool] = True,
                     do_calibrate: Optional[bool] = True,
                     do_average: Optional[bool] = True,
                     do_merge: Optional[bool] = True):
    """Combines all the facettes of data reduction into one executable function that takes
    a list of epochs and different datasets to be reduced via the MATISSE pipeline, then
    calibrated, merged and averaged

    Parameters
    ----------
    data_dir: Path
    stem_dir: Path
    target_dirs: List[Path]
    do_reduce: bool, optional
    do_calibrate: bool, optional
    do_average: bool, optional
    do_merge: bool, optional
    """
    for target_dir in target_dirs:
        start_time = time.time()
        array = "ATs" if "ATs" in target_dir else "UTs"
        cprint(f"Starting data improvement of {target_dir}!", "r")
        cprint("----------------------------------------------------------------------",
              "lg")
        if do_reduce:
            reduction_pipeline(data_dir, stem_dir, target_dir, array)
        if do_calibrate:
            calibration_pipeline(data_dir, stem_dir, target_dir,)
        if do_average:
            averaging_pipeline(data_dir, stem_dir, target_dir)
        if do_merge:
            merging_pipeline(data_dir, stem_dir, target_dir)
        cprint("----------------------------------------------------------------------",
              "lg")
        cprint(f"Reduction, calibration, merging and averaging complete in "\
               f"{datetime.timedelta(seconds=(time.time()-start_time))} hh:mm:ss",
              "r")
    print("----------------------------------------------------------------------")
    print("All done!")


if __name__ == "__main__":
    data_dir = "~/Data/hd163296/"
    stem_dir = ""
    target_lst = ["20190323"]
    target_dirs = [os.path.join("ATs", folder) for folder in target_lst]
    matadrs_pipeline(data_dir, stem_dir, target_dirs)

