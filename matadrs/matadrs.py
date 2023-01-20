from pathlib import Path
from typing import List, Union, Optional

from .libs.reduction import reduce, calibrate, average, merge
from .libs.utils.tools import cprint, print_execution_time


@print_execution_time
def matadrs_pipeline(raw_dirs: Union[List[Path], [Path]],
                     product_dir: Path,
                     do_reduce: Optional[bool] = True,
                     do_calibrate: Optional[bool] = True,
                     do_average: Optional[bool] = True,
                     do_merge: Optional[bool] = True) -> None:
    """Combines all the facettes of data reduction into one executable function that takes
    a single or a list of epochs to be reduced via the MATISSE pipeline, then calibrates,
    merges and averages them, in succession

    Parameters
    ----------
    observation_dirs: List[Path] | Path
    product_dir: Path
        The directory to contain the reduced, calibrated, averaged, and merged files
    do_reduce: bool, optional
    do_calibrate: bool, optional
    do_average: bool, optional
    do_merge: bool, optional

    Notes
    -----
    WARNING: All files in a given folder from any previous reduction will be REMOVED.
    Subdirectories for the individual steps ('reduced', 'calib', 'averaged', 'final')
    will be AUTOMATICALLY created in the 'product_dir'
    """
    if isinstance(raw_dirs, list):
        if not all(list(map(lambda x: isinstance(x, Path), raw_dirs))):
            raw_dir = list(map(Path, raw_dirs))
        if not all(list(map(lambda x: x.exists(), raw_dirs))):
            raise IOError("Paths that have been input do not exists!")
    elif isinstance(raw_dirs, Path):
        if raw_dirs.exists():
            raw_dirs = [raw_dirs]
    else:
        raise IOError("Neither valid Lists of Paths nor valid Path has been input!")

    for raw_dir in raw_dirs:
        cprint(f"Starting data reduction of '{raw_dir}'...", "cy")
        cprint(f"{'':-^50}", "lg")
        if do_reduce:
            reduce(raw_dir)
        if do_calibrate:
            calibrate(data_dir, stem_dir, target_dir,)
        if do_average:
            average(data_dir, stem_dir, target_dir)
        if do_merge:
            merge(data_dir, stem_dir, target_dir)
        cprint(f"{'':-^50}", "lg")
    cprint(f"{'':-^50}", "lg")
    cprint("All done!", "cy")


if __name__ == "__main__":
    data_dir = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/matisse/GTO/hd142666")
    observation_dirs = ["UTs/20190514", "UTs/20220420", "UTs/20220422"]
    matadrs_pipeline(list(map(lambda x: data_dir / x, observation_dirs)),
                     do_reduce=True, do_calibrate=False, do_average=False, do_merge=False)
