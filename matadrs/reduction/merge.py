"""  """
from typing import List
from pathlib import Path

from .avg_oifits import oifits_patchwork
from ..utils.plot import Plotter
from ..utils.readout import ReadoutFits
from ..utils.tools import cprint, split_fits

__all__ = ["get_output_file_path", "merge_averaged_files",
           "merge_non_averaged_files", "merge_folders", "merge"]

HEADER_TO_REMOVE = [{'key': 'HIERARCH ESO INS BCD1 ID', 'value': ' '},
                    {'key': 'HIERARCH ESO INS BCD2 ID', 'value': ' '},
                    {'key': 'HIERARCH ESO INS BCD1 NAME', 'value': ' '},
                    {'key': 'HIERARCH ESO INS BCD2 NAME', 'value': ' '}]


OI_TYPES = [['flux'], ['visamp'], ['visphi'], ['vis2'], ['t3']]


# TODO: Remove the folder's input and just get the data from the fits-file?
def get_output_file_path(fits_file: Path, output_dir: Path, pipeline: bool = False) -> Path:
    """Makes the output file's Path from the information from the input (.fits)-file
    and the directory it is contained in.

    Parameters
    ----------
    fits_file : pathlib.Path
    output_dir : pathlib.Path
    pipeline : bool, optional

    Returns
    -------
    output_file : pathlib.Path
    """
    readout = ReadoutFits(fits_file)
    tpl_start_cal, detector, tpl_start_sci = str(fits_file.parent.name).split(".")[2:-1]
    if pipeline:
        return output_dir / f"{readout.name}_{tpl_start_cal}:{tpl_start_sci}_{detector}_FINAL_TARGET_INT_pipeline.fits"
    else:
        return output_dir / f"{readout.name}_{tpl_start_cal}:{tpl_start_sci}_{detector}_FINAL_TARGET_INT.fits"


def merge_files(directory: Path, ):
    coherent_flux, incoherent_flux = [directory / flux for directory in directories]
    _, incoherent_vis = [directory / vis for directory in directories]
    coherent_bcd_vis, _ = [directory / bcd for directory in directories]
    coherent_bcd_pip_vis, incoherent_bcd_pip_vis = [directory / bcd_pip for
                                                    directory in directories]
    out_file = get_output_file_path(coherent_flux, output_dir)

    flux_chopped, vis_chopped = map(lambda x: x.replace("INT", "INT_CHOPPED"),
                                   [flux, vis])
    bcd_chopped, bcd_pip_chopped = map(lambda x: x.replace("INT", "INT_CHOPPED"),
                                      [bcd, bcd_pip])

    # NOTE: The files in the 'files_to_merge' list correspond to the 'OI_TYPE' list. Thus
    # one can determine what is merged
    if "HAWAII" in str(directories[0]):
        files_to_merge = [incoherent_flux, coherent_flux,
                          coherent_bcd_vis, incoherent_vis,
                          coherent_bcd_vis]
    else:
        files_to_merge = [incoherent_flux, coherent_flux,
                          coherent_bcd_pip_vis,
                          incoherent_bcd_pip_vis,
                          coherent_bcd_pip_vis]


def merge_averaged_files(directories: List[Path], output_dir: Path) -> None:
    """Merges the averaged files of visibility-, flux- and bcd-calibration.

    Parameters
    ----------
    directories : list of pathlib.Path
        List of the Paths to the coherent- and incoherent directory.
    output_dir : pathlib.Path
    """
    # TODO: If pipeline does not have a BCD-output use calibBCD2
    # TODO: Make this into one function that takes care of both chopped and unchopped case
    flux, vis = "TARGET_AVG_FLUX_INT.fits", "TARGET_AVG_VIS_INT.fits"
    bcd, bcd_pip = "TARGET_AVG_T3PHI_INT.fits", "TARGET_CAL_INT_noBCD.fits"

    try:
        oifits_patchwork(list(map(str, files_to_merge_unchopped)), str(out_file_unchopped),
                         oi_types_list=OI_TYPES, headerval=HEADER_TO_REMOVE)
    except Exception:
        pass

    try:
        oifits_patchwork(list(map(str, files_to_merge_chopped)), str(out_file_chopped),
                         oi_types_list=OI_TYPES, headerval=HEADER_TO_REMOVE)
    except Exception:
        pass


def merge_non_averaged_files(coherent_dir: Path,
                             incoherent_dir: Path, output_dir: Path) -> None:
    """Merges the non-averaged individual files of the visibility and flux
    calibration.

    Parameters
    ----------
    coherent_dir : pathlib.Path
        The coherently reduced and calibrated directories.
    incoherent_dir : pathlib.Path
        The incoherently reduced and calibrated directories.
    output_dir : pathlib.Path
    """
    output_dir = output_dir / "non_averaged"
    flux_tag, vis_tag = "TARGET_FLUX_CAL", "TARGET_CAL"
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    coherent_unchopped_vis_fits, coherent_chopped_vis_fits =\
        split_fits(coherent_dir, vis_tag)
    coherent_unchopped_flux_fits, coherent_chopped_flux_fits =\
        split_fits(coherent_dir, flux_tag)
    incoherent_unchopped_vis_fits, incoherent_chopped_vis_fits =\
        split_fits(incoherent_dir, vis_tag)
    incoherent_unchopped_flux_fits, incoherent_chopped_flux_fits =\
        split_fits(incoherent_dir, flux_tag)

    for index, (coh_unchopped_vis, inc_unchopped_vis, coh_unchopped_flux, inc_unchopped_flux)\
            in enumerate(zip(coherent_unchopped_vis_fits, incoherent_unchopped_vis_fits,
                             coherent_unchopped_flux_fits, incoherent_unchopped_flux_fits), start=1):
        out_file_unchopped = get_output_file_path(coh_unchopped_vis, output_dir)
        out_file_unchopped = out_file_unchopped.parent / f"{out_file_unchopped.stem}_00{index}.fits"

        if "lband" in str(coherent_dir):
            files_to_merge_unchopped = [inc_unchopped_flux, coh_unchopped_flux,
                                        coh_unchopped_vis, inc_unchopped_vis,
                                        coh_unchopped_vis]
        else:
            files_to_merge_unchopped = [inc_unchopped_flux, coh_unchopped_flux,
                                        coh_unchopped_vis, inc_unchopped_vis,
                                        coh_unchopped_vis]

        oifits_patchwork(list(map(str, files_to_merge_unchopped)),
                         str(out_file_unchopped), oi_types_list=OI_TYPES)
    try:
        for index, (coh_chopped_vis, inc_chopped_vis, coh_chopped_flux, inc_chopped_flux)\
                in enumerate(zip(coherent_chopped_vis_fits, incoherent_chopped_vis_fits,
                                 coherent_chopped_flux_fits, incoherent_chopped_flux_fits), start=1):
            out_file_chopped = get_output_file_path(coh_chopped_vis, output_dir)
            out_file_chopped = out_file_chopped.parent / f"{out_file_chopped.stem}_00{index}.fits"

            if "lband" in str(coherent_dir):
                files_to_merge_chopped = [inc_chopped_flux, coh_chopped_flux,
                                          coh_chopped_vis, inc_chopped_vis,
                                          coh_chopped_vis]
            else:
                files_to_merge_chopped = [inc_chopped_flux, coh_chopped_flux,
                                          coh_chopped_vis, inc_chopped_vis,
                                          coh_chopped_vis]
            oifits_patchwork(list(map(str, files_to_merge_chopped)),
                             str(out_file_chopped), oi_types_list=OI_TYPES)
    except Exception:
        pass


def merge_folders(coherent_dirs: List[Path],
                  incoherent_dirs: List[Path], output_dir: Path) -> None:
    """Merges the averaged files of the given folders as well as the non-averaged files,
    which are looked up.

    Parameters
    ----------
    coherent_dirs : list of pathlib.Path
    incoherent_dir : list of pathlib.Path
    output_dir : Path
    """
    cprint("Merging files...", "lp")
    for coherent_dir, incoherent_dir in zip(coherent_dirs, incoherent_dirs):
        cprint(f"Merging files of folder {coherent_dir.name.split('/')[~0]}...",
               "lp")
        cprint(f"{'':-^50}", "lg")
        cprint("Merging averaged files...", "g")
        merge_averaged_files([coherent_dir, incoherent_dir], output_dir)
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
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'.
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
        unwrap = True if "HAWAII" not in fits_file.name else False
        Plotter(fits_file, save_path=output_dir).add_mosaic(unwrap=unwrap).plot(save=True, error=True)
    for fits_file in (output_dir / "non_averaged").glob("*.fits"):
        plot_fits = Plotter(fits_file, save_path=(output_dir / "non_averaged"))
        plot_fits.add_uv().add_vis(corr_flux=True).add_vis2().add_cphases(unwrap=True).plot(save=True, error=True)
    cprint(f"{'':-^50}", "lg")
    cprint("Merging Done!", "lp")
    cprint(f"{'':-^50}", "lg")
