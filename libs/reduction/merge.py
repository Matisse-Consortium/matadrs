from pathlib import Path

# TODO: Find way to make this into a complete module -> More pythonic!
from plot import Plotter
from utils import cprint, oifits_patchwork


def merge_vis_and_cphases(stem_dir: Path, average_dir: Path) -> str:
    """Merges the vis and cphases files in the respective directory

    Parameters
    ----------
    stem_dir: Path
    average_dir: Path

    Returns
    -------
    out_file: str
    """
    target_name = stem_dir.split("/")[~1]
    epoch = average_dir.name.split(".")[2]
    lband = True if "HAWAII" in average_dir else False
    band = "L" if lband else "N"
    fits_files = average_dir.glob("*.fits")
    cphases_file = [directory for directory in fits_files if "t3" in directory.lower()][0]
    vis_file  = [directory for directory in fits_files if "vis" in directory.lower()][0]
    out_file = average_dir / f"{target_name}_{epoch}_{band}_TARGET_AVG_INT.fits"
    oifits_patchwork(cphases_file, vis_file, out_file)
    # plot = Plotter([out_file], lband=lband, save_path=out_file.parent)
    # plot.add_cphases().add_corr_flux().plot(save=True)
    return out_file


# FIXME: May be wrong. Check Jozsef's code as he is doing merging differently
def merge_mode_folders(coherent_dir: Path,
                       incoherent_dir: Path, outfile_dir: Path) -> None:
        """Calls Jozsef's code and does a average over the files for one band to average
        the reduced and calibrated data

        """
        outfile_dir = outfile_dir / coherent_dir.name
        if not outfile_dir.exists():
            outfile_dir.mkdir()

        lband = True if "HAWAII" in coherent_dir else False

        incoherent_fits = incoherent_dir.glob("*.fits")
        incoherent_fits.sort(key=lambda x: x[-8:])
        coherent_fits = coherent_dir.glob("*.fits")
        coherent_fits.sort(key=lambda x: x[-8:])

        for incoherent_file, coherent_file in zip(incoherent_fits, coherent_fits):
            out_file = outfile_dir / coherent_file.name
            oifits_patchwork(incoherent_file, coherent_file, out_file)
            cprint(f"Plotting {out_file.name}")
            plot = Plotter([out_file], lband=lband, save_path=outfile_dir)
            plot.add_cphases().add_corr_flux()
            if plot.flux is not None:
                plot.add_flux()
            plot.plot(save=True)


def merge(data_dir: Path, stem_dir: Path, target_dir: Path) -> None:
    """This merges two (.fits)-files together into one, i.e. the  "incoherent"
    and the "coherent" files

    Parameters
    ----------
    data_dir: Path
        The directory in which the files are that should be merged, i.e. "incoherent" and
        "coherent" files
    stem_dir: Path
    target_dir: Path
    """
    root_dir = Path(data_dir, stem_dir, "products", target_dir)
    merge_dir = root_dir / "calib"
    outfile_dir = root_dir / "merged_and_calib"

    coherent_dirs = (merge_dir / "coherent").glob("*.rb")
    incoherent_dirs = [Path(str(directory)).replace("coherent", "incoherent")\
                       for directory in coherent_dirs]

    for coherent_dir, incoherent_dir in zip(coherent_dirs, incoherent_dirs):
        cprint("Merging (.fits)-files of folder {coherent_dir.name.split('/')[~0]}", "lp")
        cprint("------------------------------------------------------------", "lg")
        merge_mode_folders(coherent_dir, incoherent_dir, outfile_dir)
        print("------------------------------------------------------------")
        print("Done!")
        print("------------------------------------------------------------")
    cprint("Merging Done!", "lp")


if __name__ == "__main__":
    data_dir = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"
    stem_dir, target_dir = "matisse/GTO/hd163296/", "ATs/20190323"
    merge(data_dir, stem_dir, target_dir)
