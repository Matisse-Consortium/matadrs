import os
import time
import datetime

from pathlib import Path
from typing import List, Optional

from reduction import reduction_pipeline
from calibration import calibration_pipeline
from average import averaging_pipeline
from merge import merging_pipeline
from utils import cprint

# TODO: Add functionality that clears all the paths before it write again, as to overwrite them


def matadrs_pipeline(data_dir: Path, stem_dir: Path, raw_dirs: List[Path]) -> None:
    """Combines all the facettes of data reduction into one executable function that takes
    a list of epochs and different datasets to be reduced via the MATISSE pipeline, then
    calibrated, merged and averaged

    Parameters
    ----------
    raw_dirs: List[Path]
        The respective RAW directories

    See Also
    --------
    reduction_pipeline()
    calibration_pipeline()
    merging_pipeline()
    averaging_pipeline()
    """
    for target_dir in target_dirs:
        start_time = time.time()
        array = "ATs" if "ATs" in target_dir else "UTs"
        cprint(f"Starting data improvement of {target_dir}!", "r")
        cprint("----------------------------------------------------------------------",
              "lg")
        # reduction_pipeline(data_dir, stem_dir, target_dir, array)
        calibration_pipeline(data_dir, stem_dir, target_dir,)
        # averaging_pipeline(data_dir, stem_dir, target_dir)
        # merging_pipeline(data_dir, stem_dir, target_dir)
        cprint("----------------------------------------------------------------------",
              "lg")
        cprint(f"Reduction, calibration, merging and averaging complete in "\
               f"{datetime.timedelta(seconds=(time.time()-start_time))} hh:mm:ss",
              "r")
    print("----------------------------------------------------------------------")
    print("All done!")


if __name__ == "__main__":
    data_dir = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"
    stem_dir = "matisse/GTO/hd142666/"
    target_lst = ["20190514", "20220420", "20220422"]
    target_dirs = [os.path.join("UTs", folder) for folder in target_lst]
    matadrs_pipeline(data_dir, stem_dir, target_dirs)

