"""Tool for flux and correlated flux calibration of VLTI/MATISSE data.

Jozsef Varga, 2019-2022
varga@strw.leidenuniv.nl

Example usage:
--------------
>>> from fluxcal import fluxcal
>>> inputfile_sci = 'path/to/raw/science/fits/file.fits'
>>> inputfile_cal = 'path/to/raw/calibrator/fits/file.fits'
>>> outputfile = 'path/to/calibrated/outputfile.fits'
>>> cal_database_dir = 'path/to/calibrator/database/folder/'
>>> cal_database_paths = [cal_database_dir+'vBoekelDatabase.fits',cal_database_dir+'calib_spec_db_v10.fits',cal_database_dir+'calib_spec_db_v10_supplement.fits']
>>> output_fig_dir = 'path/to/figure/folder/'
>>> fluxcal(inputfile_sci, inputfile_cal, outputfile, cal_database_paths, mode='corrflux',output_fig_dir=output_fig_dir)

Arguments:
* inputfile_sci: path to the raw science oifits file.
* inputfile_cal: path to the raw calibrator oifits file.
* outputfile: path of the output calibrated file
* cal_database_paths: list of paths to the calibrator databases, e.g., [caldb1_path, caldb2_path]
* mode (optional):
  'flux': calibrates total flux (incoherently processed oifits file expected)
          results written in the OI_FLUX table (FLUXDATA column)
  'corrflux': calibrates correlated flux (coherently processed oifits file expected)
              results written in the OI_VIS table (VISAMP column)
  'both': calibrates both total and correlated fluxes
* output_fig_dir (optional): If it is a valid path, the script will make a plot of the calibrator model spectrum there,
                             Default: '' (= no figure made)
"""
# NOTE: Load vanBoekeldatabase: DONE
# NOTE: Save calibrator spectrum in a plot: DONE
# TODO: Airmass correction: testing
# TODO: Calculate uncertainties
#       partly implemented
#       Caveats: uncertainty in the calibrator spectrum is not taken into account
#                uncertainty in calibrator diameter is not taken into account
# TODO: Treat if the cal database cannot be opened
# NOTE: Treat if there is no matching source in the database: DONE
# FIXME: Some dec values are NULL in vBoekeldatabase

########################################################################
import os
import math
import subprocess
from pathlib import Path
from shutil import copyfile
from typing import Any, List, Optional, Tuple, Union
from warnings import warn

import astropy.units as u
import numpy as np
from numpy._typing import ArrayLike
import scipy.stats
from astropy.convolution import Gaussian1DKernel, Box1DKernel, convolve
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astroquery.simbad import Simbad
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from numpy.polynomial.polynomial import polyval
from scipy.special import j1
from scipy.interpolate import interp1d


# NOTE: Try to import skycalc_cli or will install it.
# Executable should then be detected on system.
try:
    import skycalc_cli
except ImportError:
    warn("Skycalc_cli not found (is needed for proper function of this script).")
    user_input = input("Do you want to install skycalc_cli (y/n): ")
    install = True if user_input.lower() in ["yes", "no"] else False
    if install:
        subprocess.check_call(['pip', 'install', 'skycalc_cli'])


def get_spectrum_caldb(cal_database_path: Path,
                       cal_name: str, out_lst: List,
                       ra: Optional[float] = np.nan,
                       dec: Optional[float] = np.nan, 
                       match_radius: Optional[float] = 20.0,
                       band: Optional[str] = "L") -> bool:
    """Gets the spectrum calibrator.

    Parameters
    ----------
    cal_database_path : pathlib.Path or str
        Path to the calibrator database.
    cal_name : str
        Name of the calibrator.
    out_lst : list of tuples
        List of tuples (cal_name_db, diam_cal, diam_err_cal, wav_cal, spectrum_cal, caldb_file, min_sep_arcsec, ra_cal_db, dec_cal_db)
    ra : float
        Right ascension of the calibrator in degrees.
    dec : float
        Declination of the calibrator.
    match_radius : float
        Match radius in degrees.
    band : str

    Returns
    -------
    match : bool
        Whether the calibrator was found.
    """
    c_cal = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    caldb = fits.open(cal_database_path)
    caldb_file = os.path.basename(cal_database_path)
    if 'vBoekelDatabase' in caldb_file:
        if 'fitsold' in caldb_file:
            cal_name_lst = caldb[8].data['NAME']
            cal_ra_lst = caldb[8].data['RAEPP']
            cal_dec_lst = caldb[8].data['DECEPP']
        else:
            cal_name_lst = caldb['SOURCES'].data['NAME']
            cal_ra_lst = caldb['SOURCES'].data['RAEPP']
            cal_dec_lst = caldb['SOURCES'].data['DECEPP']
    elif 'calib_spec_db' in caldb_file:
        cal_name_lst = caldb['SOURCES'].data['name']
        cal_ra_lst = caldb['SOURCES'].data['ra']
        cal_dec_lst = caldb['SOURCES'].data['dec']
    else:
        cal_name_lst = caldb['SOURCES'].data['NAME']
        cal_ra_lst = caldb['SOURCES'].data['RAEPP']
        cal_dec_lst = caldb['SOURCES'].data['DECEPP']

    c_lst = SkyCoord(cal_ra_lst * u.deg, cal_dec_lst * u.deg, frame='icrs')

    # NOTE: Search for the calibrator in the calibrator database
    sep = c_cal.separation(c_lst)
    min_sep_idx = np.nanargmin(sep)
    min_sep = sep[min_sep_idx]
    if (min_sep < match_radius*u.deg/3600.0):
        # NOTE: Match
        print('Calibrator found in the database '+caldb_file+': '\
                + cal_name_lst[min_sep_idx] + ', separation: %.2f arcsec' % (3600.0*min_sep.value))

        # NOTE: Get calibrator diameter
        if 'vBoekelDatabase' in caldb_file:
            offset = 9
            if 'fitsold' in caldb_file:
                diam_cal =  1000.0*caldb[-2].data['DIAMETER'][min_sep_idx] # mas
                diam_err_cal = 1000.0*caldb[-2].data['DIAMETER_ERR'][min_sep_idx] # mas
            else:
                diam_cal =  1000.0*caldb['DIAMETERS'].data['DIAMETER'][min_sep_idx] # mas
                diam_err_cal = 1000.0*caldb['DIAMETERS'].data['DIAMETER_ERR'][min_sep_idx] # mas
        elif 'calib_spec_db' in caldb_file:
            offset = 2
            if 'L' in band:
                diam_cal = caldb['SOURCES'].data['UDDL_est'][min_sep_idx] # mas
                diam_err_cal = caldb['SOURCES'].data['e_diam_est'][min_sep_idx] # mas
            if 'N' in band:
                diam_cal = caldb['SOURCES'].data['UDDN_est'][min_sep_idx] # mas
                diam_err_cal = caldb['SOURCES'].data['e_diam_est'][min_sep_idx] # mas
            if math.isnan(diam_cal):
                diam_cal = caldb['SOURCES'].data['diam_midi'][min_sep_idx] # mas
                diam_err_cal = caldb['SOURCES'].data['e_diam_midi'][min_sep_idx] #mas
            if math.isnan(diam_cal):
                diam_cal = caldb['SOURCES'].data['diam_cohen'][min_sep_idx] # mas
                diam_err_cal = caldb['SOURCES'].data['e_diam_cohen'][min_sep_idx] # mas
            if math.isnan(diam_cal):
                diam_cal = caldb['SOURCES'].data['UDD_meas'][min_sep_idx] # mas
                diam_err_cal = caldb['SOURCES'].data['e_diam_meas'][min_sep_idx] # mas
            if math.isnan(diam_cal):
                diam_cal = caldb['SOURCES'].data['diam_gaia'][min_sep_idx] # mas
                diam_err_cal = diam_cal*0.1 # mas

        # NOTE: Extract calibrator model spectrum.
        wav_cal = caldb[min_sep_idx+offset].data['WAVELENGTH'] # m
        spectrum_cal = caldb[min_sep_idx+offset].data['FLUX']
        out_lst += [caldb[min_sep_idx+offset].header['NAME'], diam_cal, diam_err_cal,
            wav_cal, spectrum_cal, caldb_file, 3600.0*min_sep.value,
            cal_ra_lst[min_sep_idx], cal_dec_lst[min_sep_idx]]
        caldb.close()
        return True
    else:
        print('Calibrator not found in '+caldb_file)
        print('Closest match: ' + cal_name_lst[min_sep_idx] + ', separation: %.2f arcsec' % (3600.0*min_sep.value))
        caldb.close()
        return False


def plot_spec(cal_name: str, cal_database_path: Path,
              fig_dir: Optional[Path] ='.',
              ra: Optional[float] = np.nan, 
              dec: Optional[float] = np.nan,
              match_radius: Optional[float] = 15.0, 
              wl_lim: Optional[Tuple[float, float]] = (np.nan, np.nan), 
              xlog: Optional[bool] = False,
              ylog: Optional[bool] = False) -> None:
    """Plots the spectrum of the calibrator.

    Parameters
    ----------
    cal_name : str
        Name of the calibrator.
    cal_database_path : pathlib.Path or str
        Path to the calibrator database.
    fig_dir : pathlib.Path or str, optional
        Path to the figure directory.
    ra : float, optional
        Right ascension of the calibrator.
    dec : float, optional
        Declination of the calibrator.
    match_radius : float, optional
        Match radius in degrees.
    wl_lim : tuple of floats, optional
        Lower and upper wavelength limits.
    xlog : bool, optional
        If True, x-axis will be logarithmic.
    ylog : bool, optional
        If True, y-axis will be logarithmic.
    """
    if math.isnan(ra):
        Simbad.reset_votable_fields()
        Simbad.remove_votable_fields('coordinates')
        Simbad.add_votable_fields('ra(d)', 'dec(d)')
        res =  Simbad.query_object(cal_name)
        ra = res['RA_d'][0]
        dec = res['DEC_d'][0]
    out_lst =[]
    match = get_spectrum_caldb(cal_database_path, cal_name, out_lst,
                               ra=ra, dec=dec, match_radius=15.0, band='L')
    if match:
        cal_name_db, diam_cal, diam_err_cal, wav_cal, spectrum_cal,\
                caldb_file, min_sep_arcsec, ra_cal_db, dec_cal_db = out_lst
        fig, ((axf)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
        plt.plot(wav_cal*1e6, spectrum_cal, '-k',lw=2.0)
        plt.ylabel('Flux (Jy)')
        plt.xlabel(r'$\lambda$ ($\mu$m)')
        plt.title(cal_name_db + r', diam = $%.2f \pm %.2f$ mas' % (diam_cal, diam_err_cal))
        if not math.isnan(wl_lim[0]):
            axf.set_xlim(wl_lim)
        else:
            wl_lim=(np.min(wav_cal*1e6),np.max(wav_cal*1e6))
        wl_idx = np.logical_and(wav_cal*1e6 > wl_lim[0], wav_cal*1e6 < wl_lim[1])
        if xlog:
            axf.set_xscale('log')
        if ylog:
            axf.set_yscale('log')
        else:
            axf.set_ylim((0.0,1.1*np.max(spectrum_cal[wl_idx])))
        fig.savefig(fig_dir + '/' + cal_name_db.replace(' ','_')\
                + '_%.1d-%.1d_um.png' % (wl_lim[0], wl_lim[1]), dpi=200)


def plot_skycalc_output(fpath_in: Path, figpath_out: Path,
                        airmass: float, pwv: float) -> None:
    """Plots the skycalc output.

    Parameters
    ----------
    fpath_in : pathlib.Path
        Path to the skycalc output file.
    figpath_out : pathlib.Path
        Path to the figure directory.
    airmass : float
        Airmass.
    pwv : float
        Percipitative water vapor.
    """
    hdul = fits.open(fpath_in)
    wl_um = hdul[1].data['lam']*1e-3
    trans = hdul[1].data['trans']
    fig, ((axs)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
    plt.plot(wl_um , trans,lw=1.5)
    plt.ylabel('Transmission')
    plt.xlabel(r'$\lambda$ ($\mu$m)')
    plt.title('Airmass = %.3f, pwv = %.2f mm' % (airmass, pwv))
    axs.set_ylim([0.0, 1.0])
    fig.savefig(figpath_out, dpi=200)
    hdul.close()


def fluxcal(inputfile_sci: Path, inputfile_cal: Path,
            outputfile: Path, cal_database_paths: List[Path], 
            mode: Optional[str] = "flux",
            output_fig_dir: Optional[Path] = '',
            match_radius: Optional[float] = 25.0,
            do_airmass_correction: Optional[bool] = False,
            calc_spectrum_offset: Optional[bool] = False) -> int:
    """The flux calibration.

    Parameters
    -----------
    inputfile_sci : pathlib.Path
        Path to the input science (.fits)-file.
    inputfile_cal : pathlib.Path
        Path to the input calibrator (.fits)-file.
    outputfile : pathlib.Path
        Path to the output file.
    cal_database_paths : list of pathlib.Path
        Path to the calibrator database.
    mode : str
        Either "flux", "corrflux" or "both".
    output_fig_dir : pathlib.Path, optional
        Path to the output figure directory.
    match_radius : float, optional
        Match radius in arcsec.
    do_airmass_correction : bool, optional
        Whether to do airmass correction.
    calc_spectrum_offset : bool, optional
        Whether to calculate the spectrum offset.

    Returns
    -------
    error_code : int
    """
    # NOTE: Create the output oifits file
    copyfile(inputfile_sci, outputfile)
    outhdul = fits.open(outputfile, mode='update')

    # NOTE: Read in the input oifits files
    try:
        inhdul_sci = fits.open(inputfile_sci)
    except FileNotFoundError:
        print('Target reduced data file not found.')
        outhdul.close()
        os.remove(outputfile)
        return 3
    try:
        inhdul_cal = fits.open(inputfile_cal)
    except FileNotFoundError:
        print('Calibrator reduced data file not found.')
        inhdul_sci.close()
        outhdul.close()
        os.remove(outputfile)
        return 2

    airmass_sci = (inhdul_sci[0].header['HIERARCH ESO ISS AIRM START'] +\
            inhdul_sci[0].header['HIERARCH ESO ISS AIRM END'])/2.0
    pwv_sci = (inhdul_sci[0].header['HIERARCH ESO ISS AMBI IWV30D START']\
            + inhdul_sci[0].header['HIERARCH ESO ISS AMBI IWV30D END'])/2.0
    sci_name = inhdul_sci['OI_TARGET'].data['TARGET'][0]

    # NOTE: Extract calibrator information
    cal_name = inhdul_cal['OI_TARGET'].data['TARGET'][0]
    ra_cal = inhdul_cal['OI_TARGET'].data['RAEP0'][0]
    dec_cal = inhdul_cal['OI_TARGET'].data['DECEP0'][0]
    airmass_cal = (inhdul_cal[0].header['HIERARCH ESO ISS AIRM START']\
            + inhdul_cal[0].header['HIERARCH ESO ISS AIRM END'])/2.0
    pwv_cal = (inhdul_cal[0].header['HIERARCH ESO ISS AMBI IWV30D START']\
            + inhdul_cal[0].header['HIERARCH ESO ISS AMBI IWV30D END'])/2.0
    seeing_cal = (inhdul_cal[0].header['HIERARCH ESO ISS AMBI FWHM START']\
            + inhdul_cal[0].header['HIERARCH ESO ISS AMBI FWHM END'])/2.0
    tau0_cal = (inhdul_cal[0].header['HIERARCH ESO ISS AMBI TAU0 START']\
            + inhdul_cal[0].header['HIERARCH ESO ISS AMBI TAU0 END'])/2.0
    tpl_start_cal = inhdul_cal[0].header['HIERARCH ESO TPL START']
    band = inhdul_cal[0].header['HIERARCH ESO DET CHIP TYPE'] #'IR-LM' or 'IR-N'
    print('SCI: %s, airmass = %.2f CAL: %s, RA = %.6f Dec = %.6f, airmass = %.2f'\
            % (sci_name, airmass_sci, cal_name, ra_cal, dec_cal, airmass_cal))
    print('CAL TPL START: ' + tpl_start_cal)

    wav_cal = inhdul_cal['OI_WAVELENGTH'].data['EFF_WAVE']  # m
    wav_sci = inhdul_sci['OI_WAVELENGTH'].data['EFF_WAVE']  # m

    if do_airmass_correction:
        print('Do airmass correction.')
        outputdir = os.path.dirname(outputfile)+'/skycalc/'
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        tag_sci = os.path.splitext(os.path.basename(inputfile_sci))[0]

        wmin = np.min(wav_sci)*1e9 # [nm]
        wmax = np.max(wav_sci)*1e9
        margin = 0.1*(wmax-wmin)
        wmin = wmin-margin
        wmax = wmax+margin
        dlambda = get_dlambda(inhdul_sci)

        fname_in_sci = outputdir + '/skycalc_input_sci_' + tag_sci + '.txt'
        fname_skycalc_out_sci = outputdir + '/skycalc_output_sci_' + tag_sci + '.fits'
        create_skycalc_inputfile(fname_in_sci, airmass_sci, pwv_sci, wmin, wmax, wdelta=dlambda)
        print('Start SkyCalc (SCI).')
        subprocess.check_call(["skycalc_cli", '-i', fname_in_sci, '-o', fname_skycalc_out_sci])
        figpath_out = outputdir+'/skycalc_output_sci_'+tag_sci+'.png'
        plot_skycalc_output(fname_skycalc_out_sci, figpath_out, airmass_sci, pwv_sci)

        wmin = np.min(wav_cal)*1e9
        wmax = np.max(wav_cal)*1e9
        margin = 0.1*(wmax-wmin)
        wmin = wmin-margin
        wmax = wmax+margin
        dlambda = get_dlambda(inhdul_cal)

        tag_cal = os.path.splitext(os.path.basename(inputfile_cal))[0]
        fname_in_cal = outputdir + '/skycalc_input_cal_' + tag_cal + '.txt'
        fname_skycalc_out_cal = outputdir + '/skycalc_output_cal_' + tag_cal + '.fits'
        create_skycalc_inputfile(fname_in_cal, airmass_cal, pwv_cal, wmin, wmax, wdelta=dlambda)
        print('Start SkyCalc (CAL).')
        subprocess.check_call(["skycalc_cli", '-i', fname_in_cal, '-o', fname_skycalc_out_cal])
        figpath_out = outputdir + '/skycalc_output_cal_' + tag_cal + '.png'
        plot_skycalc_output(fname_skycalc_out_cal, figpath_out, airmass_cal,pwv_cal)

    # NOTE: Open the calibrator database which includes the spectra of calibrators
    match = False
    out_lst = []
    i = 0

    # NOTE: For cal_database_path in cal_database_paths:
    while (match == False and i < len(cal_database_paths)):
        match = get_spectrum_caldb(cal_database_paths[i], cal_name, out_lst,
                                   ra=ra_cal, dec=dec_cal, match_radius=match_radius, band=band)
        i = i+1

    if match:
        cal_name_db, diam_cal, diam_err_cal, wav_cal_model,\
                spectrum_cal, caldb_file, min_sep_arcsec, ra_cal_db, dec_cal_db = out_lst
        if math.isnan(diam_cal):
            print('Calibrator diameter not found.')
            if mode != 'flux':
                inhdul_cal.close()
                inhdul_sci.close()

                # NOTE: Changes are written back to fits
                outhdul.flush()
                outhdul.close()
                os.remove(outputfile)
                return 4
        print('Diameter = %.1f +/- %.1f mas (%s)' % (diam_cal, diam_err_cal, cal_name_db))

        if 'calib_spec' in caldb_file:
            wav_cal_model = np.flip(wav_cal_model)
            spectrum_cal = np.flip(spectrum_cal) # Jy
        wl = np.concatenate((np.array([wav_cal_model[0] - (wav_cal_model[1] - wav_cal_model[0])]), wav_cal_model))
        wh = np.concatenate((wav_cal_model, np.array([wav_cal_model[-1] + (wav_cal_model[-1] - wav_cal_model[-2])])))
        wm = (wh + wl) / 2
        wav_bin_lower_cal = wm[:-1]
        wav_bin_upper_cal = wm[1:]
        d_wav_cal = wav_bin_upper_cal - wav_bin_lower_cal
    else:
        print('Calibrator not found in any of the databases')
        inhdul_cal.close()
        inhdul_sci.close()

        # NOTE: Changes are written back to fits
        outhdul.flush()
        outhdul.close()
        os.remove(outputfile)
        return 1

    if np.min(wav_cal) < 0.0:
        print('ERROR (fluxcal): Wavelength grid in oifits file invalid.')
        inhdul_cal.close()
        inhdul_sci.close()

        # NOTE: Changes are written back to fits
        outhdul.flush()
        outhdul.close()
        os.remove(outputfile)
        return 5

    if 'L' in band:
        wav_cal = np.flip(wav_cal)
        wav_sci = np.flip(wav_sci)



    # NOTE: Smooth and resample calibrator spectrum to match the spectral
    # resolution and wavelength grid of the data
    print('Resample the calibrator spectrum to match the data.')
    dl_coeffs = get_dl_coeffs(inhdul_cal)
    spectral_binning = get_spectral_binning(inhdul_cal)
    if not math.isnan(spectral_binning) and np.all(np.isfinite(dl_coeffs)):
        # NOTE: New way (2022 April)
        # restrict calibrator spectrum wavelength range
        idx = np.logical_and(wav_cal_model > 0.99*np.nanmin(wav_cal), wav_cal_model < 1.01*np.nanmax(wav_cal))
        spectrum_cal = spectrum_cal[idx]
        wav_cal_model = wav_cal_model[idx]
        kernel_width_px = 10.0
        if os.path.exists(output_fig_dir):
            figpath_out = output_fig_dir + '/' + 'calibrator_' + band\
                    + '_' + cal_name.replace(' ', '_') + '_spectrum_resampling.png'
            make_plot = True
        else:
            figpath_out = ''
            make_plot = False
        spectrum_cal_resampled  = transform_spectrum_to_real_spectral_resolution(
            wav_cal_model*1e6, spectrum_cal, dl_coeffs, kernel_width_px,
            wav_sci*1e6, spectral_binning, make_plot=make_plot,
            figpath_out=figpath_out, ylabel='Flux (Jy)',yscale='final')
    else:
        # NOTE: Old way
        # idx=np.logical_and(wav_cal_model > np.nanmin(wav_cal) , wav_cal_model < np.nanmax(wav_cal))
        # print(len(wav_cal_model),np.nansum(idx),len(wav_cal))
        idx=np.logical_and(wav_cal_model > np.nanmin(wav_cal) , wav_cal_model < np.nanmax(wav_cal))
        if 2.0*len(wav_cal_model) < np.nansum(idx):
            # NOTE: If the sampling of the MATISSE spectrum is much sparser than the sampling of the calibrator spectrum:
            wl = np.concatenate((np.array([wav_cal_model[0] - (wav_cal_model[1] - wav_cal_model[0])]), wav_cal_model))
            wh = np.concatenate((wav_cal_model, np.array([wav_cal_model[-1] + (wav_cal_model[-1] - wav_cal_model[-2])])))
            wm = (wh+wl)/2
            wav_bin_lower = wm[:-1]
            wav_bin_upper = wm[1:]
            d_wav = wav_bin_upper - wav_bin_lower
            # print(wav_cal_model*1e6)

            # NOTE: resample model spectrum to the wavelengths of the observation
            spectrum_cal_resampled = wav_cal_model*0.0
            flux_calibrated_sci = wav_sci*0.0
            corrflux_calibrated_sci = wav_sci*0.0
            vis_cal = wav_cal_model*0.0
            # print(wav_cal_model)
            # idx=np.logical_and(wav_cal_model<13.0e-6,wav_cal_model>8.0e-6)
            # print(wav_cal_model[idx]*1e6)
            # print(wav_bin_upper_cal,wav_bin_lower_cal)
            for i in range(len(wav_cal_model)):
                # print(wav_bin_lower[i],wav_bin_upper[i])
                wu = (wav_bin_upper_cal - wav_bin_lower[i])
                wl = (wav_bin_upper[i] - wav_bin_lower_cal)
                wi = np.where(np.logical_and(wu > 0.0,wl > 0.0))
                # print(wav_cal_model[wi])
                # print(wi)
                wi = wi[0]
                # print(wi)

                # NOTE: Sum up the spectral values within the wavelength bin with weighting
                sum = 0.0
                sum = sum + (wav_bin_upper_cal[wi[0]] - wav_bin_lower[i])*spectrum_cal[wi[0]]
                sum = sum + (wav_bin_upper[i] - wav_bin_lower_cal[wi[-1]]) * spectrum_cal[wi[-1]]
                for j in range(1,len(wi)-1):
                    sum = sum + spectrum_cal[wi[j]]*d_wav_cal[wi[j]]
                spectrum_cal_resampled[i] = sum/d_wav[i]
                # print(i,len(wi),sum,d_wav[i],spectrum_cal_resampled[i])
        else:
            #NOTE: If the sampling of the MATISSE spectrum comparable to or denser
            # than the sampling of the calibrator spectrum do an interpolation
            f = interp1d(wav_cal_model, spectrum_cal,kind='cubic')
            # print(band)
            # print(wav_cal_model)
            # print(wav_cal)
            # print( inputfile_cal)
            spectrum_cal_resampled = f(wav_sci)

    # print(wav_cal_model[idx])
    # print(wav_cal)
    # print(spectrum_cal[idx])
    # print(spectrum_cal_resampled)

    # NOTE: Plot calibrator spectrum in the data wavelength range
    if os.path.exists(output_fig_dir):
        fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
        if 2.0*len(wav_cal) < np.nansum(idx):
            plt.plot(wav_cal_model*1e6, spectrum_cal, '-',color='grey',label='Original sp.',lw=1.0,alpha=0.66)
        else:
            plt.plot(wav_cal_model*1e6, spectrum_cal, '-o',color='grey',label='Original sp.',lw=1.0,alpha=0.66)
        plt.plot(wav_sci*1e6,spectrum_cal_resampled,'-r',label='Resampled total sp.')

    if do_airmass_correction:
        print('Calculate airmass correction factor')
        kernel_width_px = 10.0

        # NOTE: Read in the saved atmospheric transmission curves
        hdulc = fits.open(fname_skycalc_out_cal)
        wl_um_cal = hdulc[1].data['lam']*1e-3
        trans_cal = hdulc[1].data['trans']
        hduls = fits.open(fname_skycalc_out_sci)
        wl_um_sci = hduls[1].data['lam']*1e-3
        trans_sci = hduls[1].data['trans']

        # NOTE: Transform the transmission to the wavelength grid of the data
        # (incl. spectral resolution & spectral binning)
        dl_coeffs = get_dl_coeffs(inhdul_sci)
        spectral_binning = get_spectral_binning(inhdul_sci)
        figpath_out = outputdir+'/skycalc_spec_transform_sci_'+tag_sci+'.png'
        trans_sci_final = transform_spectrum_to_real_spectral_resolution(
                wl_um_sci, trans_sci, dl_coeffs, kernel_width_px,
                wav_sci*1e6, spectral_binning, make_plot=True,
                figpath_out=figpath_out,ylabel='Transmission')
        dl_coeffs = get_dl_coeffs(inhdul_cal)
        spectral_binning = get_spectral_binning(inhdul_cal)
        figpath_out = outputdir + '/skycalc_spec_transform_cal_' + tag_cal + '.png'
        trans_cal_final = transform_spectrum_to_real_spectral_resolution(
                wl_um_cal, trans_cal,dl_coeffs, kernel_width_px, wav_sci*1e6,
                spectral_binning, make_plot=True, figpath_out=figpath_out, ylabel='Transmission')

        # txtpath_out = outputdir+'/skycalc__'+tag_cal+'.png'
        # wertret
        # write_spectrum_txt(trans_cal_final,wav_cal*1e6,txtpath_out)

        # NOTE: Calculate the airmass_correction_factor
        airmass_correction_factor = trans_cal_final/trans_sci_final

        # NOTE: Plot the airmass_correction_factor
        tag_sci_final = os.path.splitext(os.path.basename(outputfile))[0]
        figa, ((axc)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
        plt.plot(wav_sci*1e6, airmass_correction_factor)
        wl_min = np.min(wl_um_cal)
        wl_max = np.max(wl_um_cal)
        wl_range = wl_max - wl_min
        LM_gap_idx = np.logical_and(wav_sci*1e6 > 4.1 ,wav_sci*1e6 < 4.6)
        idx = np.logical_and(wav_sci*1e6 > wl_min+wl_range*0.1, wav_sci*1e6 < wl_max-wl_range*0.1)
        range_idx = np.logical_and(idx, ~LM_gap_idx)
        pmin = np.nanmin(airmass_correction_factor[range_idx])  # np.nanpercentile(y_new, 10.0)
        pmax = np.nanmax(airmass_correction_factor[range_idx])  # np.nanpercentile(y_new, 95.0)
        axc.set_ylim([pmin, pmax])
        axc.set_ylabel('Airmass correction factor')
        axc.set_xlabel('$\lambda$ ($\mu$m)')
        # ax1.set_ylim([0.0,1.0])
        figpath_out = outputdir + '/skycalc_airmass_correction_factor_' + tag_sci_final + '.png'
        figa.savefig(figpath_out, dpi=200)
        hduls.close()
        hdulc.close()
    else:
        airmass_correction_factor = wav_sci*0.0+1.0

    # NOTE: Calibrate total spectrum
    if mode == 'flux' or mode == 'both':
        print('Calibrate total spectrum.')

        # NOTE: Check if we have an 'OI_FLUX' table
        n_exp_sci = len(inhdul_sci['OI_VIS2'].data['VIS2DATA'])/6
        n_exp_cal = len(inhdul_cal['OI_VIS2'].data['VIS2DATA'])/6
        rp_list_sci=[]
        rp_list_cal=[]
        try:
            for j in range(len(inhdul_sci['OI_FLUX'].data['FLUXDATA'])):
                flux_raw_sci = inhdul_sci['OI_FLUX'].data['FLUXDATA'][j] # *np.exp(airmass_sci)
                flux_raw_cal = inhdul_cal['OI_FLUX'].data['FLUXDATA'][j] # *np.exp(airmass_cal)
                fluxerr_raw_sci = inhdul_sci['OI_FLUX'].data['FLUXERR'][j]
                fluxerr_raw_cal = inhdul_cal['OI_FLUX'].data['FLUXERR'][j]
                if 'L' in band:
                    flux_raw_sci = np.flip(flux_raw_sci)
                    flux_raw_cal = np.flip(flux_raw_cal)
                    fluxerr_raw_sci = np.flip(fluxerr_raw_sci)
                    fluxerr_raw_cal = np.flip(fluxerr_raw_cal)

                if do_airmass_correction:
                    if calc_spectrum_offset:
                        # NOTE: Calculate correlation between the raw spectrum,
                        # and the atmospheric transmission spectrum
                        # shift_max = int(0.025*len(trans_cal_final))
                        shift_max = int(7.0*spectral_binning)
                        rp_list_cal.append(calc_corr_offset(trans_cal_final,flux_raw_cal,shift_max))
                        rp_list_sci.append(calc_corr_offset(trans_sci_final,flux_raw_sci,shift_max))

                flux_calibrated_sci = flux_raw_sci/flux_raw_cal*spectrum_cal_resampled\
                        * n_exp_cal/n_exp_sci*airmass_correction_factor
                if 'L' in band:
                    flux_calibrated_sci = np.flip(flux_calibrated_sci)
                fluxerr_calibrated_sci = np.abs(flux_raw_sci/flux_raw_cal*n_exp_cal/n_exp_sci)* \
                    np.sqrt((fluxerr_raw_sci*n_exp_cal/n_exp_sci/flux_raw_sci)**2\
                    + (fluxerr_raw_cal/flux_raw_cal)**2) * spectrum_cal_resampled
                if 'L' in band:
                    fluxerr_calibrated_sci = np.flip(fluxerr_calibrated_sci)
                outhdul['OI_FLUX'].data['FLUXDATA'][j] = flux_calibrated_sci
                outhdul['OI_FLUX'].data['FLUXERR'][j] = fluxerr_calibrated_sci
            outhdul['OI_FLUX'].header['TUNIT5'] = 'Jy'
            outhdul['OI_FLUX'].header['TUNIT6'] = 'Jy'

            if do_airmass_correction:
                if calc_spectrum_offset:
                    # NOTE: Plot the correlation
                    plot_corr_offset(rp_list_sci,np.arange(-shift_max, +shift_max), 
                                     outputdir + '/skycalc_correlation_flux_' + tag_sci + '.png')
                    plot_corr_offset(rp_list_cal,np.arange(-shift_max, +shift_max),
                                     outputdir + '/skycalc_correlation_flux_' + tag_cal + '.png')

        except KeyError as e:
            print('No OI_FLUX table found.')

    # NOTE: Calibrate correlated spectrum
    if mode == 'corrflux' or mode == 'both':
        print('Calibrate correlated spectra.')
        rp_list_sci = []
        rp_list_cal = []
        for j in range(len(inhdul_sci['OI_VIS'].data['VISAMP'])):
            sta_index_sci = inhdul_sci['OI_VIS'].data['STA_INDEX'][j]
            corrflux_raw_sci = inhdul_sci['OI_VIS'].data['VISAMP'][j] # *np.exp(airmass_sci)
            corrfluxerr_raw_sci = inhdul_sci['OI_VIS'].data['VISAMPERR'][j]
            if 'L' in band:
                corrflux_raw_sci = np.flip(corrflux_raw_sci)
                corrfluxerr_raw_sci = np.flip(corrfluxerr_raw_sci)

            # NOTE: Find calibrator data with matching station configuration
            sta_indices_cal = inhdul_cal['OI_VIS'].data['STA_INDEX']
            for i in range(len(sta_indices_cal)):
                if ((sta_index_sci[0] == sta_indices_cal[i][0]) and (sta_index_sci[1] == sta_indices_cal[i][1])) \
                or ((sta_index_sci[0] == sta_indices_cal[i][1]) and (sta_index_sci[1] == sta_indices_cal[i][0])):
                    idx_cal = i
                    break

            corrflux_raw_cal = inhdul_cal['OI_VIS'].data['VISAMP'][idx_cal] # *np.exp(airmass_cal)
            corrfluxerr_raw_cal = inhdul_cal['OI_VIS'].data['VISAMPERR'][idx_cal]
            if 'L' in band:
                corrflux_raw_cal = np.flip(corrflux_raw_cal)
                corrfluxerr_raw_cal = np.flip(corrfluxerr_raw_cal)
            uu = inhdul_cal['OI_VIS'].data['UCOORD'][idx_cal]
            vv = inhdul_cal['OI_VIS'].data['VCOORD'][idx_cal]
            B_p = np.sqrt(uu**2 + vv**2)

            diam_cal_rad = diam_cal/1000.0/3600.0*math.pi/180.0
            spatial_frequency = B_p/wav_cal
            # NOTE: Visibilities of the calibrator (uniform disk model)
            vis_cal = 2*j1(math.pi*diam_cal_rad*spatial_frequency) / (math.pi*diam_cal_rad*spatial_frequency)
            # plt.figure()
            # plt.plot(wav_cal, vis_cal, '-b')
            # plt.show()
            if do_airmass_correction:
                if calc_spectrum_offset:
                    # NOTE: Calculate correlation between the raw spectrum,
                    # and the atmospheric transmission spectrum
                    # shift_max = int(0.025*len(trans_cal_final))
                    shift_max = int(7.0*spectral_binning)
                    if len(trans_cal_final) == len(corrflux_raw_cal):
                        rp_list_cal.append(calc_corr_offset(trans_cal_final, corrflux_raw_cal, shift_max))
                        rp_list_sci.append(calc_corr_offset(trans_sci_final, corrflux_raw_sci, shift_max))

            corrflux_calibrated_sci = corrflux_raw_sci/corrflux_raw_cal*vis_cal\
                    * spectrum_cal_resampled*airmass_correction_factor
            if 'L' in band:
                corrflux_calibrated_sci = np.flip(corrflux_calibrated_sci)
            # NOTE: Uncertainty in calibrator diameter is not taken into account
            corrfluxerr_calibrated_sci = np.abs(corrflux_raw_sci/corrflux_raw_cal)* \
                np.sqrt((corrfluxerr_raw_sci/corrflux_raw_sci)**2 + (corrfluxerr_raw_cal/corrflux_raw_cal)**2)* \
                vis_cal*spectrum_cal_resampled
            if 'L' in band:
                corrfluxerr_calibrated_sci = np.flip(corrfluxerr_calibrated_sci)
            if os.path.exists(output_fig_dir):
                plt.plot(wav_cal*1e6, vis_cal*spectrum_cal_resampled,
                         '--', label='Resampl. corr., B_p = %.2f m'%B_p)
            # plt.figure()
            # plt.plot(wav_cal, vis_cal, '-b')
            # plt.show()
            outhdul['OI_VIS'].data['VISAMP'][j] = corrflux_calibrated_sci
            outhdul['OI_VIS'].data['VISAMPERR'][j] = corrfluxerr_calibrated_sci
        outhdul['OI_VIS'].header['TUNIT5'] = 'Jy'
        outhdul['OI_VIS'].header['TUNIT6'] = 'Jy'
        if do_airmass_correction:
            # NOTE: Plot the correlation
            if calc_spectrum_offset:
                plot_corr_offset(rp_list_sci, np.arange(-shift_max, +shift_max),
                                 outputdir + '/skycalc_correlation_corrflux_' + tag_sci + '.png')
                plot_corr_offset(rp_list_cal, np.arange(-shift_max, +shift_max),
                                 outputdir + '/skycalc_correlation_corrflux_' + tag_cal + '.png')

    outhdul[0].header['HIERARCH ESO PRO CAL NAME'] = cal_name
    outhdul[0].header['HIERARCH ESO PRO CAL RA'] = (ra_cal, '[deg]')
    outhdul[0].header['HIERARCH ESO PRO CAL DEC'] = (dec_cal, '[deg]')
    outhdul[0].header['HIERARCH ESO PRO CAL AIRM'] = airmass_cal
    outhdul[0].header['HIERARCH ESO PRO CAL IWV'] = pwv_cal
    outhdul[0].header['HIERARCH ESO PRO CAL FWHM'] = (seeing_cal, '[arcsec]')
    outhdul[0].header['HIERARCH ESO PRO CAL TAU0'] = (tau0_cal, 'Coherence time [s]')
    outhdul[0].header['HIERARCH ESO PRO CAL TPL START'] = tpl_start_cal
    outhdul[0].header['HIERARCH ESO PRO CAL DB NAME'] = (cal_name_db, 'Name of calibrator in cal database.')
    outhdul[0].header['HIERARCH ESO PRO CAL DB DBNAME'] = (caldb_file, 'Name of cal database')
    outhdul[0].header['HIERARCH ESO PRO CAL DB RA'] = (ra_cal_db, '[deg]')
    outhdul[0].header['HIERARCH ESO PRO CAL DB DEC'] = (dec_cal_db , '[deg]')
    outhdul[0].header['HIERARCH ESO PRO CAL DB DIAM'] = (diam_cal, 'Calibrator diameter [mas]')
    outhdul[0].header['HIERARCH ESO PRO CAL DB ERRDIAM'] = (diam_err_cal, 'Error in calibrator diameter [mas]')
    outhdul[0].header['HIERARCH ESO PRO CAL DB SEP'] = (min_sep_arcsec, 'Separation (input coord - calDB coord) [arcsec]')
    outhdul[0].header['HIERARCH ESO PRO CATG'] = 'TARGET_FLUXCAL_INT'

    if os.path.exists(output_fig_dir):
        ax1.set_ylabel('Flux (Jy)')
        ax1.set_xlabel(r'$\lambda$ ($\mu$m)')
        ax1.set_title(cal_name)
        ax1.set_xlim([np.min(wav_cal)*1e6, np.max(wav_cal)*1e6])
        ax1.set_ylim([0.7*np.min(spectrum_cal_resampled), 1.1*np.max(spectrum_cal_resampled)])
        ax1.legend(loc='best', fontsize=8, fancybox=True, framealpha=0.5)
        fig.savefig(output_fig_dir + '/' + 'calibrator_' + band\
                + '_' + cal_name.replace(' ','_') + '_spectrum.png', dpi=200)
        # plt.show()

    inhdul_cal.close()
    inhdul_sci.close()

    # NOTE: Changes are written back to fits
    outhdul.flush()  
    outhdul.close()
    return 0


def update_corrflux_from_vis2(
        inputfile: Path, outputfile: Path, flux_Jy: int) -> None:
    """Updates the correlated fluxes from the squared visibilities.

    Parameters
    ----------
    inputfile: pathlib.Path
        Input (.fits)-file.
    outputfile: pathlib.Path
        Output (.fits)-file.
    flux_Jy: int
        Flux in Jansky.
    """
    # NOTE: Create the output oifits file
    print(os.path.basename(inputfile))
    copyfile(inputfile, outputfile)
    inhdul = fits.open(inputfile)
    outhdul = fits.open(outputfile, mode='update')
    for j in range(len(outhdul['OI_VIS'].data['VISAMP'])):
        vis = np.sqrt(inhdul['OI_VIS2'].data['VIS2DATA'][j])
        viserr = np.abs(0.5*inhdul['OI_VIS2'].data['VIS2ERR'][j]/np.sqrt(inhdul['OI_VIS2'].data['VIS2DATA'][j]))
        corrflux = vis*flux_Jy
        corrfluxerr = viserr*flux_Jy
        outhdul['OI_VIS'].data['VISAMP'][j] = corrflux
        outhdul['OI_VIS'].data['VISAMPERR'][j] = corrfluxerr

    # NOTE: Changes are written back to fits
    outhdul.flush()
    outhdul.close()

def update_vis2_from_corrflux(
        inputfile: Path, outputfile: Path, total_flux) -> None:
    """Updates the squared visibilities from the correlated fluxes.

    Parameters
    ----------
    inputfile: pathlib.Path
    outputfile: pathlib.Path
    total_flux: int
    """
    # NOTE: Create the output oifits file
    print(os.path.basename(inputfile))
    copyfile(inputfile, outputfile)
    inhdul = fits.open(inputfile)
    outhdul = fits.open(outputfile, mode='update')
    for j in range(len(outhdul['OI_VIS2'].data['VIS2DATA'])):
        corrflux = inhdul['OI_VIS'].data['VISAMP'][j]
        corrfluxerr = inhdul['OI_VIS'].data['VISAMPERR'][j]
        vis = corrflux/total_flux
        viserr = corrfluxerr/total_flux
        outhdul['OI_VIS2'].data['VIS2DATA'][j] = vis**2
        outhdul['OI_VIS2'].data['VIS2ERR'][j] = 2.0*vis*viserr

    # NOTE: Changes are written back to fits
    outhdul.flush()
    outhdul.close()


def create_skycalc_inputfile(
        fname: Path, airmass: float, pwv: float, wmin: float, wmax: float,
        wgrid_mode: Optional[str] = 'fixed_wavelength_step',
        wdelta: Optional[float] = 0.1, wres: Optional[float] = 20000.0,
        lsf_type: Optional[str] = 'none', lsf_gauss_fwhm: Optional[float] = 5.0) -> None:
    """Creates a skycalc input file.

    Parameters
    ----------
    fname : pathlib.Path
    airmass : float
    pwv : float
    wmin : float
    wmax : float
    wgrid_mode : str, optional
        Either 'fixed_wavelength_step' or 'fixed_spectral_resolution'.
    wdelta : float, optional
    wres : float, optional
    lsf_type : str, optional
       Either 'none' or 'Gaussian' or 'Boxcar'.
    lsf_gauss_fwhm : float, optional
    """
    pwv_allowed_values = [0.05, 0.1, 0.25, 0.5, 1.0, 1.5,
                          2.5, 3.5, 5.0, 7.5, 10.0, 20.0, 30.0]
    # print(pwv)
    ofile = open(fname,'w')
    mystr = 'airmass         :  ' + '%f\n' % airmass
    mystr += 'pwv_mode        :  pwv \n'
    mystr += 'season          :  0 \n'
    mystr += 'time            :  0 \n'
    idx = find_nearest_idx(pwv_allowed_values, pwv)
    # print(idx)
    mystr += 'pwv             :  ' + '%f\n' % pwv_allowed_values[idx]
    mystr += 'msolflux        :  130.0\n'
    mystr += 'incl_moon       :  Y\n'
    mystr += 'moon_sun_sep    :  90.0\n'
    mystr += 'moon_target_sep :  45.0\n'
    mystr += 'moon_alt        :  45.0\n'
    mystr += 'moon_earth_dist :  1.0\n'
    mystr += 'incl_starlight  :  Y\n'
    mystr += 'incl_zodiacal   :  Y\n'
    mystr += 'ecl_lon         :  135.0\n'
    mystr += 'ecl_lat         :  90.0\n'
    mystr += 'incl_loweratm   :  Y\n'
    mystr += 'incl_upperatm   :  Y\n'
    mystr += 'incl_airglow    :  Y\n'
    mystr += 'incl_therm      :  N\n'
    mystr += 'therm_t1        :  0.0\n'
    mystr += 'therm_e1        :  0.0\n'
    mystr += 'therm_t2        :  0.0\n'
    mystr += 'therm_e2        :  0.0\n'
    mystr += 'therm_t3        :  0.0\n'
    mystr += 'therm_e3        :  0.0\n'
    mystr += 'vacair          :  vac\n'
    mystr += 'wmin            :  ' + '%f\n' % wmin
    mystr += 'wmax            :  ' + '%f\n' % wmax
    mystr += 'wgrid_mode      :  ' + '%s\n' % wgrid_mode
    mystr += 'wdelta          :  ' + '%f\n' % wdelta
    mystr += 'wres            :  ' + '%f\n' % wres
    mystr += 'lsf_type        :  ' + '%s\n' % lsf_type
    mystr += 'lsf_gauss_fwhm  :  ' + '%f\n' % lsf_gauss_fwhm
    mystr += 'lsf_boxcar_fwhm :  5.0\n'
    mystr += 'observatory     :  paranal'
    ofile.write(mystr)
    ofile.close()


def get_dlambda(inhdul) -> float:
    """Gets the dlambda in nm."""
    header = inhdul[0].header
    dl = np.nan
    if 'AQUARIUS' in header['HIERARCH ESO DET CHIP NAME']:
        # band = 'N'
        dispname = header['HIERARCH ESO INS DIN NAME']
        if 'LOW' in dispname:
            dl = 30.0
        if 'HIGH' in dispname:
            dl = 3.0
    if 'HAWAII' in header['HIERARCH ESO DET CHIP NAME']:
        # band = 'LM'
        dispname = header['HIERARCH ESO INS DIL NAME']
        if 'LOW' in dispname:
            dl = 8.0
        if 'MED' in dispname:
            dl = 0.6
        if 'HIGH' in dispname:
            if '+' in dispname:
                np.nan #??????????????
            else:
                dl = np.nan #??????????????
    return dl


def get_spectral_binning(inhdul: fits.BinTableHDU) -> int:
    """Gets the spectral binning."""
    header = inhdul[0].header
    spectral_binning = np.nan
    for i in range(1, 17):
        hdrkey = 'HIERARCH ESO PRO REC1 PARAM' + '%d' % i + ' NAME'
        if hdrkey in header:
            if 'spectralBinning' in header[hdrkey]:
                hdrkey2 = 'HIERARCH ESO PRO REC1 PARAM' + '%d' % i + ' VALUE'
                spectral_binning = float(header[hdrkey2])
                # print(spectral_binning)
    return spectral_binning


def get_dl_coeffs(inhdul: fits.BinTableHDU) -> List[float]:
    """Gets the dl coefficients."""
    header = inhdul[0].header
    dl_coeffs = [ np.nan]*4
    if 'AQUARIUS' in header['HIERARCH ESO DET CHIP NAME']:
        # band = 'N'
        dispname = header['HIERARCH ESO INS DIN NAME']
        if 'LOW' in dispname:
            dl_coeffs = [ 0.10600484, 0.01502548, 0.00294806, -0.00021434]
        if 'HIGH' in dispname:
            dl_coeffs = [-8.02282965e-05, 3.83260266e-03,
                         7.60090459e-05, -4.30753848e-07]
    if 'HAWAII' in header['HIERARCH ESO DET CHIP NAME']:
        # band = 'LM'
        dispname = header['HIERARCH ESO INS DIL NAME']
        if 'LOW' in dispname:
            dl_coeffs = [ 0.09200542, -0.03281159,  0.02166703, -0.00309248]
        if 'MED' in dispname:
            dl_coeffs = [ 2.73866174e-10, 2.00286100e-03,
                         1.33829137e-06, -4.46578231e-10]
        if 'HIGH' in dispname:
            if '+' in dispname:
                dl_coeffs = [ np.nan]*4 #!!!!!!!!!!!!!!!!!!!!! needs to be implemented
            else:
                dl_coeffs = [-1.08178909e-04, 6.44793559e-04,
                             1.30502477e-04, -3.70606692e-06]
    return dl_coeffs


def transform_spectrum_to_real_spectral_resolution(
        wl_orig, spec_orig, dl_coeffs,
        kernel_width_px, wl_final,
        spectral_binning: int,
        make_plot: Optional[bool] = False,
        figpath_out: Optional[Path] = 'spectrum_plots.png',
        ylabel: Optional[str] = '',
        yscale: Optional[str] = 'none') -> np.ndarray:
    """Transforms the spectrum to real spectral resolution.

    Parameters
    ----------
    wl_orig : numpy.array
        Wavelengths of the original spectrum.
    spec_orig : numpy.array
        Spectrum of the original spectrum.
    dl_coeffs : numpy.array
        Coefficients of the dispersion law.
    kernel_width_px : float
        Width of the Gaussian kernel in pixels.
    wl_final : numpy.array
        Wavelengths of the final spectrum.
    spectral_binning : int
        Number of spectral bins.
    make_plot : bool, optional
        Whether to make a plot.
    figpath_out : pathlib.Path, optional
        Path to the figure.
    ylabel : str, optional
        Y-axis label.
    yscale : str, optional
        Y-axis scale, either 'none' or 'final' or '0-1'.

    Returns
    -------
    spectrum : numpy.array
    """
    # NOTE: Make an uneven wavelength grid
    # print('make an uneven wavelength grid')
    min_wl = np.min(wl_orig)
    max_wl = np.max(wl_orig)
    # print(min_wl,max_wl)
    wl_new=[]
    wl = min_wl
    wl_new.append(wl)
    while wl < max_wl:
        wl = wl + polyval(wl,dl_coeffs)/kernel_width_px
        wl_new.append(wl)
    wl_new = np.array(wl_new)
    # print(wl_orig,wl_new)

    # NOTE: Interpolate the original spectrum to the new grid
    f_spec_new = interp1d(wl_orig, spec_orig,kind='cubic',fill_value='extrapolate')
    spec_new = f_spec_new(wl_new)

    # NOTE: Convolve with Gaussian kernel
    kernel = Gaussian1DKernel(stddev=kernel_width_px/(2.0*np.sqrt(2.0*np.log(2.0))))
    spec_new[0] = np.nanmedian(spec_new[0:int(kernel.dimension/2.0)])
    spec_new[-1] = np.nanmedian(spec_new[-1:-int(kernel.dimension/2.0)])
    spec_convolved = convolve(spec_new, kernel, boundary='extend')

    # NOTE: Interpolate the convolved spectrum to the input wavelength grid
    f_spec_new = interp1d(wl_new, spec_convolved, kind='cubic', fill_value='extrapolate')
    spec_interp = f_spec_new(wl_final)

    # NOTE: Apply spectral binning: convolve with a top-hat kernel of size spectral_binning
    if spectral_binning > 1:
        kernel = Box1DKernel(spectral_binning)
        spec_final = convolve(spec_interp, kernel, boundary='extend')
    else:
        spec_final = spec_interp

    if make_plot:
        fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(20, 10))
        plt.plot(wl_orig, spec_orig, label='1) original sp.')
        plt.plot(wl_new, spec_new, label='2) regridded sp.')
        plt.plot(wl_new, spec_convolved, label=r'3) convolved sp.'+r' ($R \approx %.0f$)'\
                % np.nanmean(wl_final/polyval(wl_final, dl_coeffs)))
        plt.plot(wl_final, spec_interp, label='4) convolved interp sp.')
        plt.plot(wl_final, spec_final, label='5) binned sp.'+' (%d px)' % spectral_binning)
        plt.ylabel(ylabel)
        plt.xlabel(r'$\lambda$ ($\mu$m)')
        ax1.set_xlim([np.nanmin(wl_final), np.nanmax(wl_final)])
        if yscale == 'final':
            ax1.set_ylim([0.95*np.nanmin(spec_final), 1.05*np.nanmax(spec_final)])
        if yscale == '0-1':
            ax1.set_ylim([0.0, 1.0])
        plt.legend()
        fig.savefig(figpath_out, dpi=200)
    return spec_final


def calc_corr_offset(spectrum1: ArrayLike,
                     spectrum2: ArrayLike, 
                     shift_max: int) -> List[float]:
    """Calculate the correlation offset.

    Parameters
    ----------
    spectrum1 : array_like
        First spectrum.
    spectrum2 : array_like
        Second spectrum.
    shift_max : int
        Maximum shift.

    Returns
    -------
    rp : list of float
        List of correlation offsets.
    """
    Ntr  = len(spectrum1)
    rp = []
    for k in range(-shift_max, +shift_max):
        if k < 0:
            # print(j,Ntr,len(trans_cal_final[0:(Ntr+k)]),len(corrflux_raw_cal[-k:]))
            rp.append(scipy.stats.pearsonr(spectrum1[0:(Ntr+k)], spectrum2[-k:])[0]) #Pearson's r
        else:
            rp.append(scipy.stats.pearsonr(spectrum1[k:], spectrum2[0:(Ntr-k)])[0]) #Pearson's r
    return rp


def plot_corr_offset(rp_list: List[float], x: List[float],
                     figpath_out: Union[str, Path]) -> None:
    """ Plot the correlation offset.

    Parameters
    ----------
    rp_list : list of float
        List of correlation offsets.
    x : list of float
        List of wavelength offsets.
    figpath_out : str or pathlib.Path
        Path to save the plot.
    """
    if len(rp_list) > 0:
        fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
        for rp in rp_list:
            plt.step(x, np.array(rp),where='mid')
        plt.ylabel('Correlation')
        plt.xlabel('Wavelength offset (px)')
        tlocs, tlabels = plt.xticks()
        tspacing = tlocs[1]-tlocs[0]
        if tspacing < 20.0:
            minor_locator = MultipleLocator(1)
        elif tspacing < 40.0:
            minor_locator = MultipleLocator(2)
        elif tspacing < 80.0:
            minor_locator = MultipleLocator(5)
        else:
            minor_locator = MultipleLocator(10)
        ax1.xaxis.set_minor_locator(minor_locator)
        plt.grid(axis='x',color='0.95',which='minor')
        plt.grid(axis='x',color='0.8',which='major')
        fig.savefig(figpath_out, dpi=200)


def find_nearest_idx(array: ArrayLike, value: Any) -> int:
    """Find the index of the nearest value in an array.

    Parameters
    ----------
    array : array_like
    value : Any

    Returns
    -------
    index : int
    """
    return (np.abs(np.asarray(array) - value)).argmin()
