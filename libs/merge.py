import os

from glob import glob
from pathlib import Path
from shutil import copyfile
from astropy.io import fits
from typing import Any, Dict, List, Union, Optional

"""Slight rewrite of Jozsef's code and folder wide application"""

# The data path to the general data
# DATA_PATH = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/"

def oifits_patchwork(incoherent_file: str, coherent_file: str,
                     outfile_path: str,
                     oi_types_list=[['vis2','visamp','visphi','t3','flux']],
                     headerval=[]) -> None:
    """Jozsef's file to merge two (.fits)-files, slightly reworked"""
    if os.path.exists(incoherent_file):
        copyfile(incoherent_file, outfile_path)
    else:
        raise RuntimeError('ERROR (oifits_patchwork): File not found: '+incoherent_file)

    outhdul  = fits.open(outfile_path, mode='update')

    n_oi_types_list = len(oi_types_list)
    for i in range(n_oi_types_list):
        oi_types = oi_types_list[i]
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
                    pass

            infile2 = coherent_file
            inhdul2 = fits.open(infile2, mode='readonly')

            outhdul['OI_VIS'].header['AMPTYP'] = inhdul2['OI_VIS'].header['AMPTYP']
            outhdul['OI_VIS'].data = inhdul2['OI_VIS'].data

            #look up visphi
            if 'visphi' in oi_types:
                visphi = inhdul2['OI_VIS'].data['VISPHI']
                visphierr = inhdul2['OI_VIS'].data['VISPHIERR']
                #match station indices
                sta_index_visamp = outhdul['OI_VIS'].data['STA_INDEX']
                sta_index_visamp = [ list(item) for item in sta_index_visamp ]
                sta_index_visphi = inhdul2['OI_VIS'].data['STA_INDEX']
                sta_index_visphi = [ list(item) for item in sta_index_visphi ]
                for k in range(len(sta_index_visamp)):
                    for l in range(len(sta_index_visphi)):
                        if ((sta_index_visamp[k] == sta_index_visphi[l]) \
                            or (sta_index_visamp[k][::-1] == sta_index_visphi[l] )):
                            outhdul['OI_VIS'].data['VISPHI'][k] = inhdul2['OI_VIS'].data['VISPHI'][l]
                            outhdul['OI_VIS'].data['VISPHIERR'][k] = inhdul2['OI_VIS'].data['VISPHIERR'][l]

    for dic in headerval:
        del outhdul[0].header[dic['key']]
        outhdul[0].header[dic['key']] = dic['value']


    outhdul.flush()  # changes are written back to original.fits
    outhdul.close()
    inhdul.close()
    inhdul2.close()

def single_merge(merge_dir: Path, bands: List[str]) -> None:
    """Makes a

    Parameters
    ----------
    merge_dir: Path
        The folder in which the calibrated data is contained
    bands: List[str]
        A list containing the bands to be averaged over
    """
    for i in bands:
        incoherent_dir = glob(os.path.join(merge_dir,
                                               f"incoherent/{i}/calib", "*.rb"))
        coherent_dir = [os.path.join(merge_dir,
                                         f"coherent/{i}/calib", os.path.basename(x))\
                            for x in incoherent_dir]

        outfile_dir = os.path.join(merge_dir, "combined", i)

        for j, k in enumerate(incoherent_dir):
            print(f"Merging {os.path.basename(k)} with "\
                  f"{os.path.basename(coherent_dir[j])}")
            print("------------------------------------------------------------")
            outfile_dir = os.path.join(merge_dir, "combined",\
                                       i, os.path.basename(k))


            if not os.path.exists(outfile_dir):
                os.makedirs(outfile_dir)

            incoherent_fits_files = glob(os.path.join(k, "*.fits"))
            incoherent_fits_files.sort(key=lambda x: x[-8:])
            coherent_fits_files = glob(os.path.join(coherent_dir[j], "*.fits"))
            coherent_fits_files.sort(key=lambda x: x[-8:])

            for l, m in enumerate(incoherent_fits_files):
                print(f"Processing {os.path.basename(m)} with "\
                      f"{os.path.basename(coherent_fits_files[l])}")
                outfile_path = os.path.join(outfile_dir, os.path.basename(m))
                oifits_patchwork(m, coherent_fits_files[l], outfile_path)

            print("------------------------------------------------------------")
            print("Done!")
            print("------------------------------------------------------------")

def merging_pipeline(merge_dir: Path,
                     both: Optional[bool] = False,
                     lband: Optional[bool] = False) -> None:
    """This merges two (.fits)-files together into one, i.e. the  "incoherent"
    and the "coherent" files

    Parameters
    ----------
    merge_dir: Path
        The directory in which the files are that should be merged, i.e. "incoherent" and
        "coherent" files
    both: bool, optional
        If both bands are to be merged, this has to be false for the "lband" option to be
        considered
    lband: bool, optional
        If "both=False" and this is "True"", then lband will be merged , if
        "both=False" and this is "False", then nband will be merged
    """
    if both:
        bands = ["lband", "nband"]
    else:
        bands = ["lband"] if lband else ["nband"]

    single_merge(merge_dir, bands)

if __name__ == "__main__":
    specific_dir = "/Users/scheuck/Documents/data/20190514"
    # specific_dir = "GTO/hd142666/PRODUCT/UTs/20190514"
    merging_pipeline(specific_dir)
    # merging_pipeline(os.path.join(DATA_PATH, specific_dir), both=True)

