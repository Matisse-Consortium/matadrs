import os

from glob import glob
from pathlib import Path
from shutil import copyfile
from typing import Optional, List

from astropy.io import fits


def cprint(message: str, c: Optional[str] = None) -> None:
    """Prints with color"""
    color_dict = {"r": ["\033[91m", "\033[00m"], "g": ["\033[92m", "\033[00m"],
                  "y": ["\033[93m", "\033[00m"], "lp": ["\033[94m", "\033[00m"],
                  "p": ["\033[95m", "\033[00m"], "cy": ["\033[96m", "\033[00m"],
                  "lg": ["\033[97m", "\033[00m"]}

    if c:
        colored_string = color_dict[c]
        colored_string.insert(1, message)
        print("".join(colored_string))
    else:
        print(message)


def check_if_target(target_dir: Path) -> bool:
    """Checks if the given path contains TAR-files

    Parameters
    ----------
    target_dir: Path
        The directory that is to be checked for 'TARGET_RAW_INT' files

    Returns
    -------
    contains_target: bool
    """
    return True if glob(os.path.join(target_dir, "TARGET_RAW_INT*")) else False


def get_path_descriptor(root_dir: Path, descriptor: Path,
                        tar_dir: Path, cal_dir: Path) -> Path:
    """Assembles the names for the new directories that will contain the
    calibrated files and returns the 'output_dir'

    Parameters
    ----------
    root_dir: Path
        The root directory of the PRODUCT
    descriptor: Path
        The desired name 'TAR-CAL' directory
    tar_dir: Path
        The target directory
    cal_dir: Path
        The calibrator directory

    Returns
    -------
    output_dir: Path
        The path for the output directory
    """
    mode_and_band = str(tar_dir.parents[1]).split("/")[-2:]
    dir_name, time_stamp_sci, detector = str(tar_dir.parent).split(".")[:-1]
    dir_name = dir_name.split('/')[-1].replace("raw", "cal")
    time_stamp_cal = str(cal_dir.parent).split('.')[-3]
    new_dir_name = '.'.join([descriptor, dir_name,
                             time_stamp_sci, detector, time_stamp_cal, "rb"])
    return root_dir / "calib" / mode_and_band[0] / new_dir_name


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
                    cprint("No 'oi_flux' has been found!", "y")

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


if __name__ == "__main__":
    print(cprint("Hello", "lp"))
