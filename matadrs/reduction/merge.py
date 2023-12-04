from collections import namedtuple
from typing import Optional, List
from pathlib import Path

import numpy as np

from .avg_oifits import oifits_patchwork
from ..utils.plot import Plotter
from ..utils.readout import ReadoutFits
from ..utils.tools import cprint, split_fits

HEADER_TO_REMOVE = [{'key': 'HIERARCH ESO INS BCD1 ID', 'value': ' '},
                    {'key': 'HIERARCH ESO INS BCD2 ID', 'value': ' '},
                    {'key': 'HIERARCH ESO INS BCD1 NAME', 'value': ' '},
                    {'key': 'HIERARCH ESO INS BCD2 NAME', 'value': ' '}]


# NOTE: The files in the 'files_to_merge' list correspond to the 'OI_TYPE' list. Thus
# one can determine what is merged
OI_TYPES = [['flux'], ['visamp'], ['visphi'], ['vis2'], ['t3']]


def get_output_file_path(fits_file: Path, output_dir: Path,
                         chopped: Optional[bool] = False) -> Path:
    """Makes the output file's Path from the information from the input (.fits)-file
    and the directory it is contained in.

    Parameters
    ----------
    fits_file : pathlib.Path
        The input file.
    output_dir : pathlib.Path
    chopped : bool, optional

    Returns
    -------
    output_file : pathlib.Path
        The output file's path.
    """
    readout = ReadoutFits(fits_file)
    tpl_start_cal, detector, tpl_start_sci = str(fits_file.parent.name).split(".")[2:-1]
    if chopped:
        return output_dir / f"{readout.name}_{tpl_start_cal}:{tpl_start_sci}_{detector}_FINAL_TARGET_INT_CHOPPED.fits"
    return output_dir / f"{readout.name}_{tpl_start_cal}:{tpl_start_sci}_{detector}_FINAL_TARGET_INT.fits"


# TODO: Rewrite this at some point so it only returns really existing files?
# Other code would then need to be changed as well
# TODO: Check for the pipeline chopped as well -> See how that is done
def get_averaged_files(directories: List[Path], chopped: Optional[bool] = False):
    """Gets the averaged files' paths."""
    flux, vis = "TARGET_AVG_FLUX_INT.fits", "TARGET_AVG_VIS_INT.fits"
    bcd, bcd_pip = "TARGET_BCD_CAL_T3PHI_INT.fits", "TARGET_CAL_INT_noBCD.fits"

    if chopped:
        flux, vis = map(lambda x: x.replace("INT", "INT_CHOPPED"), [flux, vis])
        bcd = bcd.replace("INT", "INT_CHOPPED")

    fluxes = [directory / flux for directory in directories]
    visibilities = [directory / vis for directory in directories]
    bcd_visibilities = [directory / bcd for directory in directories]
    bcd_pip_visibilities = [directory / bcd_pip for directory in directories]
    return fluxes, visibilities, bcd_visibilities, bcd_pip_visibilities


def execute_multiple_merges(fluxes: List[Path], visibilities: List[Path],
                            chopped: bool, output_dir: Path) -> None:
    """Executes a merge for multiple, non-averaged files."""
    try:
        for index, observables in enumerate(zip(*fluxes, *visibilities), start=1):
            flux, vis = observables[:2], observables[2:]
            execute_merge(output_dir, flux, vis, chopped=chopped, index=index)
    except Exception:
        pass


def prepare_multiple_merges(directories: Path, output_dir: Path):
    """Gets the non-averaged (.fits)-files and executes multiple merges with it."""
    flux_tag, vis_tag = "TARGET_FLUX_CAL", "TARGET_CAL_INT_0"
    coherent_dir, incoherent_dir = map(lambda x: Path(str(x).replace("averaged", "calib").replace("-AVG", "")),
                                       directories)
    coherent_unchopped_vis, coherent_chopped_vis=\
        split_fits(coherent_dir, vis_tag)
    coherent_unchopped_flux, coherent_chopped_flux=\
        split_fits(coherent_dir, flux_tag)
    incoherent_unchopped_vis, incoherent_chopped_vis=\
        split_fits(incoherent_dir, vis_tag)
    incoherent_unchopped_flux, incoherent_chopped_flux=\
        split_fits(incoherent_dir, flux_tag)

    execute_multiple_merges([coherent_unchopped_flux, incoherent_unchopped_flux],
                            [coherent_unchopped_vis, incoherent_unchopped_vis],
                            False, output_dir)
    execute_multiple_merges([coherent_chopped_flux, incoherent_chopped_flux],
                            [coherent_chopped_vis, incoherent_chopped_vis],
                            True, output_dir)


# TODO: If pipeline does not have a BCD-output use calibBCD2
# TODO: Add phase flip for older pipelines
def execute_merge(output_dir: Path,
                  fluxes: List[Path],
                  visibilities: List[Path],
                  bcd_visibilities: Optional[List[Path]] = None,
                  bcd_pip_visibilities: Optional[List[Path]] = None,
                  chopped: Optional[bool] = False,
                  index: Optional[int] = None) -> None:
    """Prepares the lists for the 'oifits_patchwork' and executes a merge."""
    coherent_flux, incoherent_flux = fluxes
    coherent_vis, incoherent_vis = visibilities
    out_file = get_output_file_path(coherent_flux, output_dir, chopped)
    MergeFiles = namedtuple("MergeFiles", np.array(OI_TYPES.copy()).flatten())
    MergeFiles.flux, MergeFiles.visamp = incoherent_flux, coherent_flux
    MergeFiles.vis2 = incoherent_vis
    MergeFiles.visphi = MergeFiles.t3 = coherent_vis

    if index is not None:
        out_file = out_file.parent / f"{out_file.stem}_00{index}.fits"

    if "HAWAII" in str(fluxes[0]):
        if bcd_visibilities is not None:
            coherent_bcd_vis, incoherent_bcd_vis = bcd_visibilities
            MergeFiles.visphi = coherent_bcd_vis
            MergeFiles.t3 = incoherent_bcd_vis
        else:
            MergeFiles.t3 = incoherent_vis
    else:
        if bcd_visibilities is not None:
            coherent_bcd_vis, incoherent_bcd_vis = bcd_visibilities
            MergeFiles.vis2 = incoherent_bcd_vis
            MergeFiles.visphi = MergeFiles.t3 = coherent_bcd_vis
        if bcd_pip_visibilities is not None:
            coherent_pip_vis, incoherent_pip_vis = bcd_pip_visibilities
            MergeFiles.vis2 = incoherent_pip_vis
            MergeFiles.visphi = MergeFiles.t3 = coherent_pip_vis
    files_to_merge = [MergeFiles.flux, MergeFiles.visamp,
                      MergeFiles.visphi, MergeFiles.vis2, MergeFiles.t3]
    oifits_patchwork(list(map(str, files_to_merge)), str(out_file),
                     oi_types_list=OI_TYPES, headerval=HEADER_TO_REMOVE)


# TODO: Find way to remove the try statements
# TODO: Make these errors more precise at some future point and log them
def merge_averaged_files(coherent_dir: Path,
                         incoherent_dir: Path, output_dir: Path) -> None:
    """Merges the averaged files of visibility-, flux- and bcd-calibration.

    Parameters
    ----------
    directories : list of pathlib.Path
        List of the Paths to the coherent- and incoherent directory.
    output_dir : pathlib.Path
    """
    try:
        averaged_non_chopped_files = get_averaged_files([coherent_dir, incoherent_dir])
        execute_merge(output_dir, *averaged_non_chopped_files)
    except Exception as e:
        print("Error merging non-chopped files, skipping file!\n", e)
        raise e
        pass

    try:
        averaged_chopped_files = get_averaged_files(
            [coherent_dir, incoherent_dir], chopped=True)
        if averaged_chopped_files[0][0].exists():
            execute_merge(output_dir, *averaged_chopped_files, chopped=True)
    except Exception as e:
        print("Error merging chopped files, skipping file!\n", e)
        pass


# TODO: Find way to remove the try statements
def merge_non_averaged_files(coherent_dir: Path,
                             incoherent_dir: Path, output_dir: Path) -> None:
    """Merges the non-averaged individual files of the visibility and flux
    calibration.

    Parameters
    ----------
    coherent_dir : pathlib.Path
        A coherently reduced and calibrated directory.
    incoherent_dir : pathlib.Path
        A incoherently reduced and calibrated directory.
    output_dir : pathlib.Path
        The output directory.
    """
    output_dir = output_dir / "non_averaged"
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    prepare_multiple_merges([coherent_dir, incoherent_dir], output_dir)


def merge_folders(coherent_dirs: List[Path],
                  incoherent_dirs: List[Path], output_dir: Path) -> None:
    """Merges the averaged files of the given folders.

    Parameters
    ----------
    coherent_dirs : list of pathlib.Path
        The coherently reduced and calibrated directories.
    incoherent_dir : list of pathlib.Path
        The incoherently reduced and calibrated directories.
    output_dir : pathlib.Path
        The output directory.
    """
    cprint("Merging files...", "lp")
    for coherent_dir, incoherent_dir in zip(coherent_dirs, incoherent_dirs):
        cprint("Merging files of folder"
               f" {coherent_dir.name.split('/')[~0]}...", "lp")
        cprint(f"{'':-^50}", "lg")
        cprint("Merging averaged files...", "g")
        merge_averaged_files(coherent_dir, incoherent_dir, output_dir)
        cprint("Merging non-averaged files...", "g")
        merge_non_averaged_files(coherent_dir, incoherent_dir, output_dir)


def merging_pipeline(averaged_dir: Path) -> None:
    """This merges two (.fits)-files together into one, i.e. the  "incoherent"
    and the "coherent" files.

    Parameters
    ----------
    averaged_dir : pathlib.Path
        The directory containing the averaged products.
    mode : str, optional
        The mode in which the reduction is to be executed.
        Either "coherent", "incoherent" or "both".
    """
    coherent_dirs = list((averaged_dir / "averaged" / "coherent").glob("*.rb"))
    incoherent_dirs = [Path(str(directory).replace("coherent", "incoherent"))\
                       for directory in coherent_dirs]

    output_dir = averaged_dir / "final"
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    merge_folders(coherent_dirs, incoherent_dirs, output_dir)
    cprint("Plotting files...", "g")
    for fits_file in output_dir.glob("*.fits"):
        unwrap = "HAWAII" not in fits_file.name
        plot_fits = Plotter(fits_file, save_path=output_dir)
        plot_fits.add_mosaic(unwrap=unwrap).plot(save=True, error=True)
    for fits_file in (output_dir / "non_averaged").glob("*.fits"):
        plot_fits = Plotter(fits_file, save_path=output_dir / "non_averaged")
        plot_fits.add_uv().add_vis(corr_flux=True)
        plot_fits.add_vis2().add_cphases(unwrap=True)
        plot_fits.plot(save=True, error=True)
    cprint(f"{'':-^50}", "lg")
    cprint("Merging Done!", "lp")
    cprint(f"{'':-^50}", "lg")
