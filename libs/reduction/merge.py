from pathlib import Path

# TODO: Find way to make this into a complete module -> More pythonic!
from plot import Plotter
from readout import ReadoutFits
from utils import cprint, get_fits
from avg_oifits import oifits_patchwork


HEADER_TO_REMOVE = [{'key':'HIERARCH ESO INS BCD1 ID','value':' '},
                    {'key':'HIERARCH ESO INS BCD2 ID','value':' '},
                    {'key':'HIERARCH ESO INS BCD1 NAME','value':' '},
                    {'key':'HIERARCH ESO INS BCD2 NAME','value':' '}]

OI_TYPES = [["flux"], ["visamp"], ["visphi"], ["vis2"], ["t3"]]


def get_output_file_path(coherent_dir: Path, target_name: str, output_dir: Path) -> Path:
    """Makes the output file's Path from the information in the folder it is in"""
    tpl_start_cal, detector, tpl_start_sci = str(coherent_dir.name).split(".")[2:-1]
    return output_dir / f"{target_name}_{tpl_start_cal}:{tpl_start_sci}_FINAL_TARGET_INT.fits"


def merge_averaged_files(coherent_dir: Path,
                         incoherent_dir: Path, output_dir: Path) -> None:
    """Merges the averaged files of visibility-, flux- and bcd-calibration

    Parameters
    ----------
    coherent_dir: Path
    incoherent_dir: Path
    output_dir: Path
    """
    coherent_flux = str(coherent_dir / "TARGET_AVGFLUXCAL_INT.fits")
    coherent_bcd_vis = str(coherent_dir / "TARGET_AVG_T3PHI_INT.fits")
    incoherent_flux = str(incoherent_dir / "TARGET_AVGFLUXCAL_INT.fits")
    incoherent_vis = str(incoherent_dir / "TARGET_AVGCAL_INT.fits")
    incoherent_bcd_vis = str(incoherent_dir / "TARGET_AVG_T3PHI_INT.fits")
    out_file = get_output_file_path(coherent_dir,
                                    ReadoutFits(Path(coherent_flux)).target_name,
                                    output_dir)

    if "lband" in str(coherent_dir):
        files_to_merge = [incoherent_flux, coherent_flux,
                          coherent_bcd_vis, incoherent_vis, incoherent_bcd_vis]
    else:
        files_to_merge = [incoherent_flux, coherent_flux,
                          coherent_bcd_vis, incoherent_vis, coherent_bcd_vis]
    oi_types_list = [["flux"], ["visamp"], ["visphi"], ["vis2"], ["t3"]]
    oifits_patchwork(files_to_merge, str(out_file),
                     oi_types_list=OI_TYPES, headerval=HEADER_TO_REMOVE)
    # TODO: Implement handling of chopped files
    # TODO: Add plotting here


def merge_non_averaged_files(coherent_dir: Path,
                             incoherent_dir: Path, output_dir: Path) -> None:
    """Merges the non-averaged individual files of the visibility and flux calibration

    Parameters
    ----------
    coherent_dir: Path
    incoherent_dir: Path
    output_dir: Path
    """
    output_dir = output_dir / "non_averaged"
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    coherent_unchopped_vis_fits,\
            coherent_chopped_vis_fits = get_fits(coherent_dir, "TARGET_CAL")
    coherent_unchopped_flux_fits,\
            coherent_chopped_flux_fits = get_fits(coherent_dir, "TARGET_FLUXCAL")
    incoherent_unchopped_vis_fits,\
            incoherent_chopped_vis_fits = get_fits(incoherent_dir, "TARGET_CAL")
    incoherent_unchopped_flux_fits,\
            incoherent_chopped_flux_fits = get_fits(incoherent_dir, "TARGET_FLUXCAL")

    for index, (coh_unchopped_vis, inc_unchopped_vis, coh_unchopped_flux, inc_unchopped_flux)\
            in enumerate(zip(coherent_unchopped_vis_fits, incoherent_unchopped_vis_fits,
                             coherent_unchopped_flux_fits, incoherent_unchopped_flux_fits)):
        out_file = get_output_file_path(coherent_dir,
                                        ReadoutFits(Path(coh_unchopped_vis)).target_name,
                                        output_dir)
        out_file = out_file.parent / f"{out_file.stem}_00{index}.fits"

        if "lband" in str(coherent_dir):
            files_to_merge = [inc_unchopped_flux, coh_unchopped_flux,
                              coh_unchopped_vis, inc_unchopped_vis, inc_unchopped_vis]
        else:
            files_to_merge = [inc_unchopped_flux, coh_unchopped_flux,
                              coh_unchopped_vis, inc_unchopped_vis, coh_unchopped_vis]

        oifits_patchwork(files_to_merge, str(out_file), oi_types_list=OI_TYPES)
        # TODO: Implement handling of chopped files
        # TODO: Add plotting here


def merge_folders(coherent_dir: Path, incoherent_dir: Path, output_dir: Path):
    """Merges the averaged files of the given folders as well as the non-averaged files,
    which are looked up

    Parameters
    ----------
    coherent_dir: Path
    incoherent_dir: Path
    output_dir: Path
    """
    cprint(f"Merging files of folder {coherent_dir.name.split('/')[~0]}...", "lp")
    cprint(f"{'':-^50}", "lg")
    cprint(f"Merging averaged files...", "g")
    cprint(f"{'':-^50}", "lg")
    merge_averaged_files(coherent_dir, incoherent_dir, output_dir)
    cprint(f"{'':-^50}", "lg")
    cprint(f"Merging non-averaged files...", "g")
    merge_non_averaged_files(coherent_dir, incoherent_dir, output_dir)
    cprint(f"{'':-^50}", "lg")


def merge(root_dir: Path, stem_dir: Path, target_dir: Path) -> None:
    """This merges two (.fits)-files together into one, i.e. the  "incoherent"
    and the "coherent" files

    Parameters
    ----------
    root_dir: Path
    stem_dir: Path
    target_dir: Path
    """
    root_dir = Path(root_dir, stem_dir, "products", target_dir)
    coherent_dirs = (root_dir / "bcd_and_averaged" / "coherent").glob("*.rb")
    breakpoint()
    incoherent_dirs = [Path(str(directory).replace("coherent", "incoherent"))\
                       for directory in coherent_dirs]

    output_dir = root_dir / "final"
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    cprint(f"Merging files...", "lp")
    cprint(f"{'':-^50}", "lg")
    # FIXME: This is not executed -> Why?
    for coherent_dir, incoherent_dir in zip(coherent_dirs, incoherent_dirs):
        merge_folders(coherent_dir, incoherent_dir, output_dir)
    cprint("Merging Done!", "lp")


if __name__ == "__main__":
    root_dir = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"
    stem_dir, target_dir = "matisse/GTO/hd163296/", "ATs/20190323"
    merge(root_dir, stem_dir, target_dir)
