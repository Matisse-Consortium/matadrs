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

OI_TYPES = [["flux"], ["visamp"], ["visphi"], ["vis2"], ["t3"]]


# TODO: Remove the folder's input and just get the data from the fits-file?
def get_output_file_path(fits_file: Path, output_dir: Path, pipeline: bool = False) -> Path:
    """Makes the output file's Path from the information from the input (.fits)-file
    and the directory it is contained in

    Parameters
    ----------
    fits_file: Path
    output_dir: Path
    pipeline: bool, optional

    Returns
    -------
    output_file: Path"""
    readout = ReadoutFits(fits_file)
    tpl_start_cal, detector, tpl_start_sci = str(fits_file.parent.name).split(".")[2:-1]
    if pipeline:
        return output_dir / f"{readout.name}_{tpl_start_cal}:{tpl_start_sci}_{detector}_FINAL_TARGET_INT_pipeline.fits"
    else:
        return output_dir / f"{readout.name}_{tpl_start_cal}:{tpl_start_sci}_{detector}_FINAL_TARGET_INT.fits"


def merge_averaged_files(directories: List[Path], output_dir: Path) -> None:
    """Merges the averaged files of visibility-, flux- and bcd-calibration

    Parameters
    ----------
    directories: List[Path]
        List of the Paths to the coherent- and incoherent directory
    output_dir: Path
    """
    flux, vis = "TARGET_AVG_FLUX_INT.fits", "TARGET_AVG_VIS_INT.fits"
    bcd, bcd_pip = "TARGET_AVG_T3PHI_INT.fits", "TARGET_CAL_INT_noBCD.fits"
    coherent_flux, incoherent_flux = [directory / flux for directory in directories]
    coherent_vis, incoherent_vis = [directory / vis for directory in directories]
    coherent_bcd_vis, incoherent_bcd_vis = [directory / bcd for directory in directories]
    coherent_bcd_pip_vis, incoherent_bcd_pip_vis = [directory / bcd_pip for directory in directories]
    out_file = get_output_file_path(coherent_flux, output_dir)
    out_file_pip = get_output_file_path(coherent_flux, output_dir, True)

    # NOTE: The files in the 'files_to_merge' list correspond to the 'OI_TYPE' list. Thus
    # one can determine what is merged in what way
    if "lband" in str(directories[0]):
        files_to_merge = [incoherent_flux, coherent_flux,
                          coherent_bcd_vis, incoherent_vis, incoherent_bcd_vis]
        files_to_merge_pip = [incoherent_flux, coherent_flux,
                              coherent_bcd_pip_vis, incoherent_vis, incoherent_bcd_pip_vis]
    else:
        files_to_merge = [incoherent_flux, coherent_flux,
                          coherent_bcd_vis, incoherent_vis, coherent_bcd_vis]
        files_to_merge_pip = [incoherent_flux, coherent_flux,
                              coherent_bcd_pip_vis, incoherent_vis, coherent_bcd_pip_vis]
    if all([fits_file.exists() for fits_file in files_to_merge]):
        oifits_patchwork(list(map(str, files_to_merge)), str(out_file),
                         oi_types_list=OI_TYPES, headerval=HEADER_TO_REMOVE)

    if all([fits_file.exists() for fits_file in files_to_merge_pip]):
        oifits_patchwork(list(map(str, files_to_merge)), str(out_file_pip),
                         oi_types_list=OI_TYPES, headerval=HEADER_TO_REMOVE)
    # TODO: Implement handling of chopped files


def merge_non_averaged_files(coherent_dir: Path,
                             incoherent_dir: Path, output_dir: Path) -> None:
    """Merges the non-averaged individual files of the visibility and flux calibration

    Parameters
    ----------
    coherent_dir: Path
        The coherently reduced and calibrated directories
    incoherent_dir: Path
        The incoherently reduced and calibrated directories
    output_dir: Path
    """
    output_dir = output_dir / "non_averaged"
    flux_tag, vis_tag = "TARGET_FLUX_CAL", "TARGET_CAL"
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    coherent_unchopped_vis_fits,\
            coherent_chopped_vis_fits = split_fits(coherent_dir, vis_tag)
    coherent_unchopped_flux_fits,\
            coherent_chopped_flux_fits = split_fits(coherent_dir, flux_tag)
    incoherent_unchopped_vis_fits,\
            incoherent_chopped_vis_fits = split_fits(incoherent_dir, vis_tag)
    incoherent_unchopped_flux_fits,\
            incoherent_chopped_flux_fits = split_fits(incoherent_dir, flux_tag)

    for index, (coh_unchopped_vis, inc_unchopped_vis, coh_unchopped_flux, inc_unchopped_flux)\
            in enumerate(zip(coherent_unchopped_vis_fits, incoherent_unchopped_vis_fits,
                             coherent_unchopped_flux_fits, incoherent_unchopped_flux_fits), start=1):
        out_file = get_output_file_path(coh_unchopped_vis, output_dir)
        out_file = out_file.parent / f"{out_file.stem}_00{index}.fits"

        if "lband" in str(coherent_dir):
            files_to_merge = [inc_unchopped_flux, coh_unchopped_flux,
                              coh_unchopped_vis, inc_unchopped_vis, inc_unchopped_vis]
        else:
            files_to_merge = [inc_unchopped_flux, coh_unchopped_flux,
                              coh_unchopped_vis, inc_unchopped_vis, coh_unchopped_vis]

        if all([fits_file.exists() for fits_file in files_to_merge]):
            oifits_patchwork(list(map(str, files_to_merge)),
                             str(out_file), oi_types_list=OI_TYPES)
        else:
            raise FileNotFoundError("Files haven't been found: 'oifits_patchwork'"
                                    " cannot be executed!")
        # TODO: Implement handling of chopped files


def merge_folders(coherent_dirs: List[Path],
                  incoherent_dirs: List[Path], output_dir: Path) -> None:
    """Merges the averaged files of the given folders as well as the non-averaged files,
    which are looked up

    Parameters
    ----------
    coherent_dirs: List[Path]
    incoherent_dir: List[Path]
    output_dir: Path
    """
    cprint(f"Merging files...", "lp")
    for coherent_dir, incoherent_dir in zip(coherent_dirs, incoherent_dirs):
        cprint(f"Merging files of folder {coherent_dir.name.split('/')[~0]}...",
               "lp")
        cprint(f"{'':-^50}", "lg")
        cprint(f"Merging averaged files...", "g")
        merge_averaged_files([coherent_dir, incoherent_dir], output_dir)
        cprint(f"Merging non-averaged files...", "g")
        merge_non_averaged_files(coherent_dir, incoherent_dir, output_dir)


def merge(averaged_dir: Path) -> None:
    """This merges two (.fits)-files together into one, i.e. the  "incoherent"
    and the "coherent" files

    Parameters
    ----------
    averaged_dir: Path
        The directory containing the averaged products
    mode: str, optional
        The mode in which the reduction is to be executed. Either 'coherent',
        'incoherent' or 'both'
    """
    coherent_dirs = list((averaged_dir / "bcd_and_averaged" / "coherent").glob("*.rb"))
    incoherent_dirs = [Path(str(directory).replace("coherent", "incoherent"))\
                       for directory in coherent_dirs]

    output_dir = averaged_dir / "final"
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    merge_folders(coherent_dirs, incoherent_dirs, output_dir)
    cprint("Plotting files...", "g")
    for fits_file in output_dir.glob("*.fits"):
        plot_fits = Plotter([fits_file], save_path=output_dir)
        plot_fits.add_cphases().add_vis().plot(save=True)
    for fits_file in (output_dir / "non_averaged").glob("*.fits"):
        plot_fits = Plotter([fits_file], save_path=(output_dir / "non_averaged"))
        plot_fits.add_cphases().add_vis().plot(save=True)
    cprint(f"{'':-^50}", "lg")
    cprint("Merging Done!", "lp")
    cprint(f"{'':-^50}", "lg")
