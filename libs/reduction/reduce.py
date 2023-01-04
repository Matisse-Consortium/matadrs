import time
import datetime
import shutil
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.vizier import Vizier
from astrop.coordinates import SkyCoord
from mat_tools import mat_autoPipeline as mp

from libs.reduction.utils import cprint
from libs.reduction.readout import ReadoutFits


DATA_DIR = Path(__file__).parent.parent.parent / "data"

JSDC_V2_CATALOG = Vizier(catalog="II/346/jsdc_v2")
SPECTRAL_BINNING = {"low": [5, 7], "high_ut": [5, 38], "high_at": [5, 98]}
JSDC_CATALOG = DATA_DIR / "jsdc_20170303.fits"
ADDITIONAL_CATALOG = DATA_DIR / "additional_calibrators_202207.fits"


# TODO: Implement this tpl start for reduction
# {'night':'2020-01-08', 'tpl_start':'2020-01-09T01:24:18', 'tel':'UTs','diL':'LOW','diN':'LOW'},  #V900Mon, SCI_V900Mon

def in_catalog(readout: ReadoutFits, radius: u.arcsec, catalog: Path):
    """Checks if calibrator is in catalog. Returns catalog if True, else None

    Parameters
    ----------
    readout: ReadoutFits
    radius: u.arcsec
    catalog_fits: Path

    Returns
    -------
    catalog: Path | None
    """
    coords_calibrator = SkyCoord(readout.ra*u.deg, readout.dec*u.deg, frame="icrs")
    table = Table().read(catalog)
    coords_catalog = SkyCoord(table["RAJ2000"], table["DEJ2000"],
                              unit=(u.hourangle, u.deg), frame="icrs")
    separation = coords_calibrator.separation(coords_catalog)
    if separation[np.nanargmin(separation)] < radius.to(u.deg):
        return catalog
    return None


def gets_catalog_match(fits_file: Path, target_name: str,
                       match_radius: u.arcsec = 20*u.arcsec):
    """Checks if the calibrator is in the 'jsdc_v2'-catalog and if not then searches the
    local calibrator databases

    Parameters
    ----------
    fits_file: Path
        The file that is to be reduced and checked if calibrator or not
    target_name: str
        The name of the calibrator
    match_radius: u.arcsec
        The radius in which to search the catalouge

    Returns
    -------
    catalog: Path | None
    """
    readout = ReadoutFits(fits_file)
    if JSDC_V2_CATALOG.query_object(readout.target_name, radius=match_radius):
        return JSDC_CATALOG
    return in_catalog(readout, radius=match_radius, catalog=ADDITIONAL_CATALOG)


def set_script_arguments(corr_flux: bool, array: str,
                         resolution: Optional[str] = "low") -> str:
    """Sets the arguments that are then passed to the 'mat_autoPipeline.py'
    script

    Parameters
    ----------
    corr_flux: bool
        This specifies if the flux is to be reduced or not
    array: str
        The array configuration that was used for the observation
    resolution: str
        This determines the spectral binning. Input can be "low" for
        low-resolution in both bands, "high_ut" for low-resolution in L-band
        and high-resolution in N-band for the UTs and the same for "high_at" for
        the ATs

    Returns
    -------
    str
        A string that contains the arguments, which are passed to the
        MATISSE-pipline
    """
    bin_L, bin_N = SPECTRAL_BINNING[resolution]
    # NOTE: Jozsef uses 3 here, but Jacob 2? What is the meaning of this -> Read up on it
    # -> Already asked, awaiting response
    compensate = "/compensate=[pb,rb,nl,if,bp,od]"
    tel = "/replaceTel=3" if array == "ATs" else "/replaceTel=0"
    coh_L  = f"/corrFlux=TRUE/useOpdMod=FALSE/coherentAlgo=2"\
            if corr_flux else ""
    coh_N = f"/corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2" if corr_flux else ""
    paramL_lst = f"{coh_L}{compensate}/spectralBinning={bin_L}"
    paramN_lst = f"{tel}{coh_N}/spectralBinning={bin_N}"
    return paramL_lst, paramN_lst


def reduce_mode_and_band(raw_dir: Path, calib_dir: Path, res_dir: Path,
                         array: str, mode: bool, band: str, tp_start: str,
                         resolution: Optional[str] = "low") -> None:
    """Reduces either the lband or the nband data for either the "coherent" or
    "incoherent" setting for a single iteration/epoch.

    Removes the old (.sof)-files to ensure proper reduction and then creates the needed
    folders in the "res_dir"-directory. After this, it starts the reduction with the
    specified settings

    Parameters
    ----------
    raw_dir: Path
        The path containing the raw observation files
    calib_dir: Path
        The path containing the calibration files
    res_dir: Path
        The path to contain to reduced data
    array: str
        The array configuration that was used for the observation. Either "AT" or "UT"
    mode: bool
        The mode in which the reduction is to be done, either "incoherent" if "False" or
        "coherent" if "True"
    band: str
        The band for which the reduction is to be done, either "lband" or "nband"
    tpl_star: str
        The starting time of target's observations
    """
    # TODO: Replace this time with process finish time
    start_time = time.time()
    mode_and_band_dir = res_dir / mode / band
    param_L, param_N = set_script_arguments(mode, array, resolution)
    skip_L = True if band == "nband" else False
    skip_N = not skip_L

    if not mode_and_band_dir.exists():
        mode_and_band_dir.mkdir()

    # TODO: Set tpl and split calib and reduction files
    mp.mat_autoPipeline(dirRaw=raw_dir, dirResult=res_dir, dirCalib=calib_dir,
                        tplstartsel=tpl_start, nbCore=6, resol='',
                        paramL=param_L, paramN=param_N, overwrite=0, maxIter=1,
                        skipL=skip_L, skipN=skip_N)

    try:
        rb_folders = res_dir.glob("Iter1/*.rb")
        for folder in rb_folders:
            if (mode_and_band_dir / folder.name).exists():
                shutil.rmtree(mode_and_band_dir / folder.name)
            shutil.move(folder, mode_and_band_dir)

        if rb_folders:
            cprint("Folders have sucessfully been moved to their directories", "g")
    # TODO: Make logger here
    except Exception:
        cprint("Moving of files to {mode_and_band_dir} failed!", "y")

    print("---------------------------------------------------------------------")
    cprint(f"Executed the {mode} reduction for the {band} in"
          f" {datetime.timedelta(seconds=(time.time()-start_time))} hh:mm:ss")
    cprint("---------------------------------------------------------------------",
          "lg")

# # TODO: Check if object is calibrator
# if fpath_CAL_lst: #if object is a calibrator
    # hdr = hdr_lst[0]

    # # NOTE: check whether calibrator is in JSDC
    # match_radius = 20.0 #arcsec
    # jsdc_match = False
    # jsdc_path_match = ''
    # target_name = hdr['HIERARCH ESO OBS TARG NAME']
    # target_ra = hdr['RA'] #J2000
    # target_dec = hdr['DEC'] #J2000
    # #print(target_ra,target_dec)
    # print('Check whether calibrator '+target_name+' is in JSDC v2:')
    # result_jsdc = jsdc_v2.query_object(target_name,catalog='II/346/jsdc_v2',radius=match_radius*u.arcsec)
    # #print(result_jsdc)
    # if result_jsdc != []:
        # if(len(result_jsdc[0]) > 0):
            # #match
            # jsdc_match = True
            # jsdc_path_match = main_JSDC_path
            # print('Calibrator found in JSDC v2: '+result_jsdc[0]['Name'][0])
            # #    +', separation: %.2f arcsec'%(3600.0*min_sep.value))
    # else:
        # #check whether the calibrator is in the supplement catalog
        # c_cal = SkyCoord(target_ra* u.deg, target_dec* u.deg, frame='icrs')
        # caldb = fits.open(alt_JSDC_path)
        # cal_name_lst = caldb[1].data['NAME']
        # cal_ra_lst = caldb[1].data['RAJ2000']
        # cal_dec_lst = caldb[1].data['DEJ2000']
        # #print(cal_ra_lst[0:10], cal_dec_lst[0:10])
        # c_lst = SkyCoord(ra=cal_ra_lst ,dec=cal_dec_lst , unit=(u.hourangle, u.deg), frame='icrs')
        # # search for the calibrator in the calibrator database
        # sep = c_cal.separation(c_lst)
        # min_sep_idx = np.nanargmin(sep)
        # min_sep = sep[min_sep_idx]
        # caldb.close()
        # if (min_sep < match_radius*u.deg/3600.0): #match_radius = 20 arcsec
            # #match
            # jsdc_match = True
            # jsdc_path_match = alt_JSDC_path
            # print('Calibrator found in the supplement catalog: '+os.path.basename(jsdc_path_match)+': '+cal_name_lst[min_sep_idx]
                # +', separation: %.2f arcsec'%(3600.0*min_sep.value))

    # if jsdc_match:
        # cats_to_remove = []
        # #check for an existing JSDC catalog in dir_calib
        # fpath_JSDC_lst,hdr_lst = check_tag(dir_calib,'JSDC_CAT',None)
        # replace_JSDC_cat = True
        # if len(fpath_JSDC_lst) > 0:
            # for fpath_JSDC,hdr in zip(fpath_JSDC_lst,hdr_lst):
                # #check if the local JSDC catalog is the same as the reference JSDC catalog
                # hdr_match = getheader(jsdc_path_match)
                # if 'DATE' in hdr_match and 'DATE' in hdr:
                    # if hdr_match['DATE'] == hdr['DATE']:
                        # #if the JSDC catalog in dir_calib is the same as the one in the list of JSDC catalogs, then keep it
                        # replace_JSDC_cat = False
                # else:
                    # cats_to_remove.append(fpath_JSDC)

        # # remove unwanted JSDC cats
        # for fpath in cats_to_remove:
            # print('Removing unwanted JSDC cat: '+fpath)
            # os.remove(fpath)

        # if replace_JSDC_cat:
            # # copy the wanted JSDC catalog to dir_calib
            # print('Copying JSDC cat '+os.path.basename(jsdc_path_match)+' to '+dir_calib)
            # shutil.copyfile(jsdc_path_match,os.path.join(dir_calib,os.path.basename(jsdc_path_match)))
    # else:
        # print('Calibrator not found in JSDC, reduced data will not contain TF2, thus visibility calibration (mat_cal_oifits) will fail.')


def reduce(root_dir: Path, stem_dir: Path,
           target_dir: Path, array: str,
           resolution: Optional[str] = "low"):
    """Runs the pipeline for the data reduction

    Parameters
    ----------
    raw_dir: Path
        The path containing the raw observation files
    array: str
        The array configuration that was used for the observation

    See Also
    --------
    single_reduction()
    """
    # TODO: Replace this time with process finish time
    overall_start_time = time.time()
    raw_dir = Path(root_dir, stem_dir, "raw", target_dir).resolve()

    if (raw_dir / "calib_files").exists():
        calib_dir = raw_dir / "calib_files"
    else:
        calib_dir = raw_dir

    res_dir = Path(root_dir, stem_dir, "products", target_dir).resolve()

    if not res_dir.exists():
        res_dir.mkdir()

    # TODO: Add in the option to not remove old reduction and make new one take an
    # addional tag after its name
    try:
        folders = res_dir.glob("*")
        for folder in folders:
            if folder.exists():
                shutil.rmtree(folder)
        cprint("Cleaned up old reduction!", "y")

    # TODO: Make logger here
    except Exception:
        cprint("Cleaning up failed!", "y")


    for mode in ["coherent", "incoherent"]:
        cprint(f"Processing {mode} reduction", "lp")
        cprint("---------------------------------------------------------------------",
              "lg")
        for band in ["lband", "nband"]:
            reduce_mode_and_band(raw_dir, calib_dir, res_dir, array, mode=mode, band=band)

    cprint(f"Executed the overall reduction in"
           f" {datetime.timedelta(seconds=(time.time()-overall_start_time))} hh:mm:ss",
          "lp")


if __name__ == "__main__":
    data_dir = Path("~", "Code/matadrs/data").expanduser()
    stem_dir, target_dir = "/hd163296/", "ATs/20190323"
    print(set_script_arguments(corr_flux=True, array="UTs"))
    # reduce(data_dir, stem_dir, target_dir, "ATs")

