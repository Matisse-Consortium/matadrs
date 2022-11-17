import os
import warnings

from glob import glob
from pathlib import Path
from shutil import copyfile
from astropy.io import fits
from typing import List, Optional

from plot import Plotter
from utils import cprint


def oifits_patchwork(incoherent_file: Path, coherent_file: Path, outfile_path: Path,
                     oi_types_list: Optional[List] = [['vis2','visamp','visphi','t3','flux']],
                     headerval: Optional[List] = []) -> None:
    """Jozsef's file to merge two (.fits)-files, slightly reworked"""
    if os.path.exists(incoherent_file):
        copyfile(incoherent_file, outfile_path)
    else:
        raise IOError(f"File not found: {incoherent_file}")

    outhdul = fits.open(outfile_path, mode='update')

    for oi_types in oi_types_list:
        inhdul = fits.open(incoherent_file, mode='readonly')

        for oi_type in oi_types:
            if oi_type == 'vis2':
                outhdul['OI_VIS2'].data = inhdul['OI_VIS2'].data
            if oi_type == 't3':
                outhdul['OI_T3'].data = inhdul['OI_T3'].data
            if oi_type == 'visamp':
                try:
                    outhdul[0].header['HIERARCH ESO PRO CAL NAME']          =          inhdul[0].header['HIERARCH ESO PRO CAL NAME']
                    outhdul[0].header['HIERARCH ESO PRO CAL RA']            =            inhdul[0].header['HIERARCH ESO PRO CAL RA']
                    outhdul[0].header['HIERARCH ESO PRO CAL DEC']           =           inhdul[0].header['HIERARCH ESO PRO CAL DEC']
                    outhdul[0].header['HIERARCH ESO PRO CAL AIRM']          =          inhdul[0].header['HIERARCH ESO PRO CAL AIRM']
                    outhdul[0].header['HIERARCH ESO PRO CAL FWHM']          =          inhdul[0].header['HIERARCH ESO PRO CAL FWHM']
                    outhdul[0].header['HIERARCH ESO PRO CAL TAU0']          =          inhdul[0].header['HIERARCH ESO PRO CAL TAU0']
                    outhdul[0].header['HIERARCH ESO PRO CAL TPL START']     =     inhdul[0].header['HIERARCH ESO PRO CAL TPL START']
                    outhdul[0].header['HIERARCH ESO PRO CAL DB NAME']       =       inhdul[0].header['HIERARCH ESO PRO CAL DB NAME']
                    outhdul[0].header['HIERARCH ESO PRO CAL DB DBNAME']     =     inhdul[0].header['HIERARCH ESO PRO CAL DB DBNAME']
                    outhdul[0].header['HIERARCH ESO PRO CAL DB RA']         =         inhdul[0].header['HIERARCH ESO PRO CAL DB RA']
                    outhdul[0].header['HIERARCH ESO PRO CAL DB DEC']        =        inhdul[0].header['HIERARCH ESO PRO CAL DB DEC']
                    outhdul[0].header['HIERARCH ESO PRO CAL DB DIAM']       =       inhdul[0].header['HIERARCH ESO PRO CAL DB DIAM']
                    outhdul[0].header['HIERARCH ESO PRO CAL DB ERRDIAM']    =    inhdul[0].header['HIERARCH ESO PRO CAL DB ERRDIAM']
                    # outhdul[0].header['HIERARCH ESO PRO CAL DB SEPARATION'] = inhdul[0].header['HIERARCH ESO PRO CAL DB SEPARATION']
                except KeyError as e:
                    print(e)

            if oi_type == 'flux':
                try:
                    outhdul['OI_FLUX'].data = inhdul['OI_FLUX'].data
                except KeyError as e:
                    warnings.warn("No 'oi_flux' has been found!")

            infile2 = coherent_file
            inhdul2 = fits.open(infile2, mode='readonly')

            outhdul['OI_VIS'].header['AMPTYP'] = inhdul2['OI_VIS'].header['AMPTYP']
            outhdul['OI_VIS'].data = inhdul2['OI_VIS'].data

            #look up visphi
            if 'visphi' in oi_types:
                # Match station indices
                sta_indicies_visamp = outhdul['OI_VIS'].data['STA_INDEX']
                sta_indicies_visamp = [list(item) for item in sta_indicies_visamp]
                sta_indicies_visphi = inhdul2['OI_VIS'].data['STA_INDEX']
                sta_indicies_visphi = [list(item) for item in sta_indicies_visphi]
                for i, sta_index_visamp in enumerate(sta_indicies_visamp):
                    for j, sta_index_visphi in enumerate(sta_indicies_visphi):
                        if ((sta_index_visamp == sta_index_visphi) \
                            or (sta_index_visamp[::-1] == sta_index_visphi)):
                            outhdul['OI_VIS'].data['VISPHI'][i] = inhdul2['OI_VIS'].data['VISPHI'][j]
                            outhdul['OI_VIS'].data['VISPHIERR'][i] = inhdul2['OI_VIS'].data['VISPHIERR'][j]
    for dic in headerval:
        del outhdul[0].header[dic['key']]
        outhdul[0].header[dic['key']] = dic['value']

    outhdul.flush()  # changes are written back to original.fits
    outhdul.close()
    inhdul.close()
    inhdul2.close()


def merge_vis_and_cphases(stem_dir: Path, average_dir: Path) -> str:
    """Merges the vis and cphases files in the respective directory"""
    target_name = stem_dir.split("/")[~1]
    epoch = os.path.basename(average_dir).split(".")[2]
    lband = True if "HAWAII" in average_dir else False
    band = "L" if lband else "N"
    fits_files = glob(os.path.join(average_dir, "*.fits"))
    cphases_file = [directory for directory in fits_files\
                    if "t3" in directory.lower()][0]
    vis_file  = [directory for directory in fits_files\
                 if "vis" in directory.lower()][0]
    out_file = os.path.join(average_dir,
                            f"{target_name}_{epoch}_{band}_TARGET_AVG_INT.fits")
    oifits_patchwork(cphases_file, vis_file, out_file)
    plot = Plotter([out_file], lband=lband, save_path=os.path.dirname(out_file))
    plot.add_cphases().add_corr_flux().plot(save=True)
    return out_file


def merging_pipeline(data_dir: Path, stem_dir: Path, target_dir: Path) -> None:
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
    merge_dir = os.path.join(root_dir, "bcd_and_averaged")
    outfile_dir = os.path.join(root_dir, "combined")
    if not os.path.exists(outfile_dir):
        os.makedirs(outfile_dir)

    incoherent_dirs = glob(os.path.join(merge_dir, "incoherent", "*.rb"))
    coherent_dirs = [dir.replace("incoherent", "coherent") for dir in incoherent_dirs]

    for coherent_dir, incoherent_dir in zip(coherent_dirs, incoherent_dirs):
        cprint("Merging incoherent and coherent files of folder"\
               f" {os.path.basename(coherent_dir).split('/')[~0]}", "lp")
        cprint("------------------------------------------------------------",
              "lg")
        coherent_file = merge_vis_and_cphases(stem_dir, coherent_dir)
        incoherent_file = merge_vis_and_cphases(stem_dir, incoherent_dir)
        out_file = os.path.basename(incoherent_file).replace("AVG", "FINAL")
        oifits_patchwork(incoherent_file, coherent_file,
                         os.path.join(outfile_dir, out_file))
        lband = True if "HAWAII" in out_file else False
        plot = Plotter([os.path.join(outfile_dir, out_file)], save_path=outfile_dir)
        plot.add_cphases().add_corr_flux().plot(save=True)
        print("------------------------------------------------------------")
        print("Done!")
        print("------------------------------------------------------------")
    cprint("Merging Done!", "lp")


if __name__ == "__main__":
    data_dir = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"
    stem_dir, target_dir = "matisse/GTO/hd163296/", "ATs/20190323"
    merging_pipeline(data_dir, stem_dir, target_dir)


