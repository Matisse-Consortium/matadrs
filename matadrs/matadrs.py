from pathlib import Path
from typing import List, Union, Optional

from .reduction import reduction_pipeline, calibration_pipeline,\
    averaging_pipeline, merging_pipeline
from .utils.tools import cprint, print_execution_time


@print_execution_time
def matadrs_pipeline(raw_dirs: Union[List[Path], Path],
                     product_dirs: Union[List[Path], Path],
                     mode: Optional[str] = "both",
                     band: Optional[str] = "both",
                     overwrite: Optional[bool] = False,
                     do_reduce: Optional[bool] = True,
                     do_calibrate: Optional[bool] = True,
                     do_average: Optional[bool] = True,
                     do_merge: Optional[bool] = True,
                     ncores: Optional[int] = 6) -> None:
    """Combines all the facettes of data reduction into one executable function
    that takes a single or a list of epochs to be reduced via the MATISSE
    pipeline, then calibrates, merges and averages them, in succession.

    Parameters
    ----------
    raw_dirs : list of Path or Path
        The directory/ies to either the raw, reduced, calibrated or averaged
        files depending on which step is selected for the pipeline, for the
        respective observation/s.
    product_dirs : list of Path or Path
        The directory/ies to contain the reduced, calibrated, averaged, and
        merged files for the respective observation/s.
    mode : str, optional
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'.
    band : str, optional
        The band in which the reduction is to be executed. Either 'lband',
        'nband' or 'both'.
    overwrite : bool, optional
        If 'True' overwrites present files from previous reduction.
    do_reduce : bool, optional
        Execute the reduction step.
    do_calibrate : bool, optional
        Execute the calibration step.
    do_average : bool, optional
        Execute the averaging step.
    do_merge : bool, optional
        Execute the merging step.
    ncores : int, optional
        The number of cores used for the data reduction.

    Notes
    -----
    WARNING : All files in a given folder from any previous reduction will be
    REMOVED.
    Subdirectories for the individual steps ("reduced", "calib", "averaged",
    "final") will be AUTOMATICALLY created in the 'product_dir'.
    """
    if isinstance(raw_dirs, list):
        if not all(list(map(lambda x: isinstance(x, Path), raw_dirs))):
            raw_dir = list(map(Path, raw_dirs))
        if not all(list(map(lambda x: x.exists(), raw_dirs))):
            raise IOError("The Paths in the raw_dirs list do not exists!")
    elif isinstance(raw_dirs, Path):
        if raw_dirs.exists():
            raw_dirs = [raw_dirs]
        else:
            raise IOError("The Path given for the raw_dirs do not exists!")
    else:
        raise IOError("Nor valid Lists of Paths nor valid Path for raw_dirs"
                      " has been input!")

    if isinstance(product_dirs, list):
        if not all(list(map(lambda x: isinstance(x, Path), product_dirs))):
            product_dirs = list(map(Path, product_dirs))
    elif isinstance(product_dirs, Path):
        product_dirs = [product_dirs]
    elif isinstance(product_dirs, str):
        product_dirs = [Path(product_dirs)]
    else:
        raise IOError("Input for product_dirs is in incorrect format!")

    for raw_dir, product_dir in zip(raw_dirs, product_dirs):
        cprint(f"Starting data reduction of '{raw_dir}'...", "cy")
        cprint(f"{'':-^50}", "lg")
        if do_reduce:
            reduction_pipeline(raw_dir, product_dir, mode, band, ncores, overwrite)
        if do_calibrate:
            calibration_pipeline(product_dir, mode, band, overwrite)
        if do_average:
            averaging_pipeline(product_dir, mode, overwrite)
        if do_merge:
            merging_pipeline(product_dir)
    cprint("Completed full data reduction!", "cy")
