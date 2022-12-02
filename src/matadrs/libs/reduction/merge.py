import os

from glob import glob
from pathlib import Path

from .plot import Plotter
from .utils import cprint, oifits_patchwork


def merge_mode_folders(coherent_dir: Path,
                       incoherent_dir: Path, outfile_dir: Path) -> None:
        """Calls Jozsef's code and does a average over the files for one band to average
        the reduced and calibrated data

        """
        outfile_dir = os.path.join(outfile_dir, os.path.basename(coherent_dir))
        if not os.path.exists(outfile_dir):
            os.makedirs(outfile_dir)

        lband = True if "HAWAII" in coherent_dir else False

        incoherent_fits = glob(os.path.join(incoherent_dir, "*.fits"))
        incoherent_fits.sort(key=lambda x: x[-8:])
        coherent_fits = glob(os.path.join(coherent_dir, "*.fits"))
        coherent_fits.sort(key=lambda x: x[-8:])

        for incoherent_file, coherent_file in zip(incoherent_fits, coherent_fits):
            out_file = os.path.join(outfile_dir, os.path.basename(coherent_file))
            oifits_patchwork(incoherent_file, coherent_file, out_file)
            cprint(f"Plotting {os.path.basename(out_file)}")
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
    root_dir = os.path.join(data_dir, stem_dir, "PRODUCTS", target_dir)
    merge_dir = os.path.join(root_dir, "calib")
    outfile_dir = os.path.join(root_dir, "merged_and_calib")

    coherent_dirs = glob(os.path.join(merge_dir, "coherent", "*.rb"))
    incoherent_dirs = [directory.replace("coherent", "incoherent")\
                       for directory in coherent_dirs]

    for coherent_dir, incoherent_dir in zip(coherent_dirs, incoherent_dirs):
        cprint("Merging (.fits)-files of folder"\
               f" {os.path.basename(coherent_dir).split('/')[~0]}", "lp")
        cprint("------------------------------------------------------------",
              "lg")
        merge_mode_folders(coherent_dir, incoherent_dir, outfile_dir)
        print("------------------------------------------------------------")
        print("Done!")
        print("------------------------------------------------------------")
    cprint("Merging Done!", "lp")


if __name__ == "__main__":
    data_dir = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"
    stem_dir, target_dir = "matisse/GTO/hd163296/", "ATs/20190323"
    merge(data_dir, stem_dir, target_dir)


