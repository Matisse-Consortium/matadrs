"""$Id: $

This file is part of the Matisse pipeline GUI series
Copyright (C) 2017- Observatoire de la Côte d'Azur

Created on Sat Mar 17 06:39:49 2018
@author: fmillour
fmillour@oca.eu

This software is a computer program whose purpose is to show oifits
files from the MATISSE instrument.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

Changelog:
2018-03-23: New functions: oi_data_select_frame, filter_oi_list, open_oi_dir,
            show_vis2_tf2_vs_time, show_oi_vs_time (jvarga)
2018-03-26: New GUI interface ready: oi_data_select_frame (jvarga)
2018-04-04: Updated GUI and extended functionality: input file/folder
            textbox, filter for target name, more bands (JHK) available
            (for e.g. AMBER data), plot with or without errorbars, plot V or
            V2 (jvarga)
2023-12-06: Updated saving of plots (automated folder creation) and documentation of functions (mScheuck)
"""
from pathlib import Path
from typing import List, Dict, Optional, Union

import matplotlib.gridspec as gridspec
import numpy as np
from astropy.time import Time
from astropy.io import fits
from matplotlib import pyplot as plt

from ..utils import robust


def smooth_tab(tab, n_smooth):
    kernel = np.ones(n_smooth)//n_smooth
    tab_smoothed = np.convolve(tab, kernel, mode='same')
    return tab_smoothed


def open_oi(oi_file: Path, band_LM: str) -> Dict:
    """Opens an oifits file and returns a dictionary with the
    relevant information."""
    try:
        hdu = fits.open(oi_file)
    except IOError:
        print(f"Unable to read fits file: {oi_file}")
        return {}

    hdr = hdu[0].header
    try:
        wl = hdu['OI_WAVELENGTH'].data['EFF_WAVE']
    except KeyError:
        return {}

    dic = {'WLEN': wl}

    dic['SEEING'] = (hdr['HIERARCH ESO ISS AMBI FWHM START']+hdr['HIERARCH ESO ISS AMBI FWHM END'])/2.
    dic['TAU0']   = (hdr['HIERARCH ESO ISS AMBI TAU0 START']+hdr['HIERARCH ESO ISS AMBI TAU0 END'])/2.
    airmass       = (hdr['HIERARCH ESO ISS AIRM START']+hdr['HIERARCH ESO ISS AIRM END'])/2.
    dic['PWV']   = airmass*(hdr['HIERARCH ESO ISS AMBI IWV START']+hdr['HIERARCH ESO ISS AMBI IWV END'])/2.

    dic['FT'] = hdr['HIERARCH ESO DEL FT SENSOR']

    target_name = hdu['OI_TARGET'].data['TARGET'][0]
    if not target_name:
        try:
            target_name = hdr['HIERARCH ESO OBS TARG NAME']
        except KeyError:
            print ("Target name not found.")
            target_name = ""
    dic['TARGET'] = target_name

    try:
        target_category = hdu['OI_TARGET'].data['CATEGORY'][0]  # "CAL" or "SCI"
    except KeyError:
        print ("Target category not found.")
        target_category = "CAL"
    dic['CATEGORY'] = target_category

    try:
        dateobs = hdr['DATE-OBS']
    except KeyError:
        dateobs = ""
    dic['DATEOBS'] = dateobs

    try:
        det_name = hdr['HIERARCH ESO DET CHIP NAME']
    except KeyError:
        print ("Detector name not found.")
        det_name = ""

    try:
        DIT = hdr["HIERARCH ESO DET SEQ1 DIT"]  # (s)
    except KeyError:
        DIT = np.nan
        print ("DIT not found")

    dic['DIT'] = DIT

    if (det_name == 'AQUARIUS'):
        band = 'N'
        try:
            flux_target=hdr['HIERARCH ESO CAT FLUX N']
        except:
            flux_target=hdr['HIERARCH ESO SEQ TARG FLUX N']
    elif (det_name == 'HAWAII-2RG'):
        band = 'LM'
        if (band_LM == 'M'):
            try:
                flux_target=hdr['HIERARCH ESO CAT FLUX M']
            except:
                flux_target=hdr['HIERARCH ESO SEQ TARG FLUX L']
        else:
            try:
                flux_target=hdr['HIERARCH ESO CAT FLUX L']
            except:
                flux_target=hdr['HIERARCH ESO SEQ TARG FLUX L']
    else:
        band = ''
    dic['BAND'] = band
    dic['TARGET_FLUX'] = flux_target

    try:
        if (det_name == 'AQUARIUS'):
            dispersion_name = hdr['HIERARCH ESO INS DIN NAME']
        elif (det_name == 'HAWAII-2RG'):
            dispersion_name = hdr['HIERARCH ESO INS DIL NAME']
    except KeyError:
        print ("Dispersion name not found.")
        dispersion_name = ""
    dic['DISP'] = dispersion_name

    try:
        BCD1 = hdr["HIERARCH ESO INS BCD1 NAME"]
        BCD2 = hdr["HIERARCH ESO INS BCD2 NAME"]
    except KeyError:
        BCD1 = ""
        BCD2 = ""
        print ("BCD NAME not found")
    dic['BCD1NAME'] = BCD1
    dic['BCD2NAME'] = BCD2
    dic['BCDCONFIG'] = BCD1+'-'+BCD2
    try:
        dic['TEL_NAME'] = hdu['OI_ARRAY'].data["TEL_NAME"]
        dic['STA_NAME'] = hdu['OI_ARRAY'].data["STA_NAME"]
        dic['STA_INDEX'] = hdu['OI_ARRAY'].data["STA_INDEX"]
    except KeyError:
        dic['TEL_NAME'] = {}
        dic['STA_NAME'] = {}
        dic['STA_INDEX'] = {}
        print ("Key in table OI_ARRAY not found")
    try:
        dic['VIS'] = {}
        dic['VIS']['VIS'] = hdu['OI_VIS'].data['VISAMP']
        dic['VIS']['VISERR'] = hdu['OI_VIS'].data['VISAMPERR']
        dic['VIS']['DPHI'] = hdu['OI_VIS'].data['VISPHI']
        dic['VIS']['DPHIERR'] = hdu['OI_VIS'].data['VISPHIERR']
        dic['VIS']['U'] = hdu['OI_VIS'].data['UCOORD']
        dic['VIS']['V'] = hdu['OI_VIS'].data['VCOORD']
        dic['VIS']['TIME'] = hdu['OI_VIS'].data['MJD']
        dic['VIS']['STA_INDEX'] = hdu['OI_VIS'].data['STA_INDEX']
    except:
        print("WARNING: No OI_VIS table!")

    try:
        dic['VIS2'] = {}
        dic['VIS2']['VIS2'] = hdu['OI_VIS2'].data['VIS2DATA']
        dic['VIS2']['VIS2ERR'] = hdu['OI_VIS2'].data['VIS2ERR']
        dic['VIS2']['U'] = hdu['OI_VIS2'].data['UCOORD']
        dic['VIS2']['V'] = hdu['OI_VIS2'].data['VCOORD']
        dic['VIS2']['TIME'] = hdu['OI_VIS2'].data['MJD']
        dic['VIS2']['STA_INDEX'] = hdu['OI_VIS2'].data['STA_INDEX']
        dic['CF2'] = {}
        n_base=6
        n_meas=np.shape(dic['VIS2']['U'])[0]
        n_lam=np.shape(dic['VIS2']['VIS2'])[1]
        n_exp=n_meas//n_base
        B_length=np.sqrt(np.power(dic['VIS2']['U'][0:6],2)+np.power(dic['VIS2']['V'][0:6],2))
        ind_ref=np.argmin(B_length)
        dic['CF2']['CF2_RATIO'] = np.zeros([n_meas,n_lam],dtype=np.float)
        dic['CF2']['CF2_RATIOERR'] = np.zeros([n_meas,n_lam],dtype=np.float)
        sigma2_CF2_rel=np.zeros([n_base,n_lam],dtype=np.float)
        sigma2_CF2_rel_ref=np.zeros(n_lam,dtype=np.float)

        for i in range(n_exp):
            for j in range(6):
                dic['CF2']['CF2_RATIO'][i*n_base+j,:] = smooth_tab(dic['VIS2']['VIS2'][i*n_base+j,:],6)/smooth_tab(dic['VIS2']['VIS2'][ind_ref+n_base*i,:],6)
                sigma2_CF2_rel=np.power(hdu['OI_VIS2'].data['VIS2ERR'][i*n_base+j,:],2)/np.power(dic['VIS2']['VIS2'][i*n_base+j,:],2)
                sigma2_CF2_rel_ref=np.power(dic['VIS2']['VIS2ERR'][ind_ref+n_base*i,:],2)/np.power(dic['VIS2']['VIS2'][ind_ref+n_base*i,:],2)
                dic['CF2']['CF2_RATIOERR'][i*n_base+j,:] = dic['CF2']['CF2_RATIO'][i*n_base+j,:]*np.sqrt(sigma2_CF2_rel+sigma2_CF2_rel_ref)
        dic['CF2']['TIME'] = dic['VIS2']['TIME']
        dic['CF2']['STA_INDEX'] = dic['VIS2']['STA_INDEX']
    except Exception:
        print("WARNING: No OI_VIS2 table!")

    try:
        dic['TF2'] = {}
        dic['TF2']['TF2'] = hdu['TF2'].data['TF2']
        dic['TF2']['TF2ERR'] = hdu['TF2'].data['TF2ERR']
        dic['TF2']['TF'] = hdu['TF2'].data['TF']
        dic['TF2']['TFERR'] = hdu['TF2'].data['TFERR']
        dic['TF2']['TIME'] = hdu['TF2'].data['MJD']
        dic['TF2']['STA_INDEX'] = hdu['TF2'].data['STA_INDEX']
    except Exception:
        print("WARNING: No OI_TF2 table!")

    try:
        dic['T3'] = {}
        dic['T3']['T3AMP'] = hdu['OI_T3'].data['T3AMP']
        dic['T3']['T3AMPERR'] = hdu['OI_T3'].data['T3AMPERR']
        dic['T3']['CLOS'] = hdu['OI_T3'].data['T3PHI']
        dic['T3']['CLOSERR'] = hdu['OI_T3'].data['T3PHIERR']
        dic['T3']['U1'] = hdu['OI_T3'].data['U1COORD']
        dic['T3']['V1'] = hdu['OI_T3'].data['V1COORD']
        dic['T3']['U2'] = hdu['OI_T3'].data['U2COORD']
        dic['T3']['V2'] = hdu['OI_T3'].data['V2COORD']
        dic['T3']['TIME'] = hdu['OI_T3'].data['MJD']
        dic['T3']['STA_INDEX'] = hdu['OI_T3'].data['STA_INDEX']
    except Exception:
        print("WARNING: No OI_T3 table!")

    try:
        dic['FLUX'] = {}
        dic['FLUX']['FLUX'] = hdu['OI_FLUX'].data['FLUXDATA']
        dic['FLUX']['FLUXERR'] = hdu['OI_FLUX'].data['FLUXERR']
        dic['FLUX']['TIME'] = hdu['OI_FLUX'].data['MJD']
        dic['FLUX']['STA_INDEX'] = hdu['OI_FLUX'].data['STA_INDEX']
    except Exception:
        print("WARNING: No OI_FLUX table!")

    hdu.close()
    return dic


def show_oi_vs_freq(dic: Dict,
                    log: Optional[bool] = False,
                    showvis: Optional[bool] = False) -> None:
    """Plots the oi vs. freq for the oifits-files.

    Parameters
    ----------
    dic : dict
        Dictionary of oifits-files.
    log : bool, optional
        If True, the plots are in log scale. The default is False.
    showvis : bool, optional
        If True, the visibility is shown. The default is False.
    """
    wl    = dic['WLEN'];
    vis2  = dic['VIS2']['VIS2'];
    vis2e = dic['VIS2']['VIS2ERR'];
    u     = dic['VIS2']['U'];
    v     = dic['VIS2']['V'];
    cp    = dic['T3']['CLOS'];
    cpe   = dic['T3']['CLOSERR'];
    u1    = dic['T3']['U1'];
    v1    = dic['T3']['V1'];
    u2    = dic['T3']['U2'];
    v2    = dic['T3']['V2'];

    plt.figure(figsize=(9, 6))
    G = gridspec.GridSpec(2, 1)

    axes_v2 = plt.subplot(G[0, :])

    # NOTE: Visibility. Plot all data first
    for i, _ in enumerate(u):
        r = np.sqrt(u[i] ** 2 + v[i] ** 2);
        freq = r / wl;
        if log:
            if showvis == True:
                axes_v2.semilogy(freq, vis2[i, :], color='lightgray')
            else:
                axes_v2.semilogy(freq, np.sqrt(vis2[i, :]), color='lightgray')
            plt.ylim([1e-4, 1.1])
        else:
            if showvis == True:
                axes_v2.plot(freq, np.sqrt(vis2[i, :]), color='lightgray')
            else:
                axes_v2.plot(freq, vis2[i, :], color='lightgray')

    # Plot valid data
    for i, _ in enumerate(u):
        r = np.sqrt(u[i] ** 2 + v[i] ** 2);
        freq = r / wl;
        test = np.logical_and(vis2[i, :] >= 0, vis2e[i, :] / vis2[i, :] < 1)
        if log:
            if showvis == True:
                axes_v2.semilogy(freq[test], np.sqrt(vis2[i, test]))
            else:
                axes_v2.semilogy(freq[test], vis2[i, test])
            plt.ylim([1e-4, 1.1])
        else:
            if showvis == True:
                axes_v2.plot(freq[test], np.sqrt(vis2[i, test]))
            else:
                axes_v2.plot(freq[test], vis2[i, test])

    plt.ylim([-0.1, 1.1])
    if showvis == True:
        plt.ylabel('V')
        axes_v2.set_title('Visibilities vs frequencies')
    else:
        plt.ylabel('V2')
        axes_v2.set_title('Squared visibilities vs frequencies')

    # NOTE: Closure phase
    # Plot all data first
    axes_cp = plt.subplot(G[1, :])
    for i, _ in enumerate(u1):
        r1 = np.sqrt(u1[i] ** 2 + v1[i] ** 2);
        r2 = np.sqrt(u2[i] ** 2 + v2[i] ** 2);
        r3 = np.sqrt((u1[i] + u2[i]) ** 2 + (v1[i] + v2[i]) ** 2);
        freq = np.maximum(np.maximum(r1, r2), r3) / wl;
        axes_cp.plot(freq, cp[i, :], color='lightgray')

    # Plot valid data only
    for i, _ in enumerate(u1):
        r1 = np.sqrt(u1[i] ** 2 + v1[i] ** 2);
        r2 = np.sqrt(u2[i] ** 2 + v2[i] ** 2);
        r3 = np.sqrt((u1[i] + u2[i]) ** 2 + (v1[i] + v2[i]) ** 2);
        freq = np.maximum(np.maximum(r1, r2), r3) / wl;
        test = np.absolute(cpe[i, :]) < 180 / np.pi / 3
        axes_cp.plot(freq[test], cp[i, test])

    plt.ylim([-200, 200])
    axes_cp.set_title('Closure phase vs frequencies')
    plt.ylabel('Closure phase')
    plt.xlabel(r'Spatial Frequency (B/$\lambda$)')
    plt.show()


def show_oi_vs_wlen(
        dic: Dict,
        key: Optional[str] = 'VIS2',
        datatype: Optional[str] = "VIS2", 
        showvis: Optional[bool] = False,
        plot_errorbars: Optional[bool] = True,
        correct_polynom: Optional[bool] = False,
        timeVertOffset: Optional[int] = 0,
        stdevRange: Optional[bool] = False,
        normRange: Optional[bool] = False,
        stdevTime: Optional[bool] = False) -> None:
    """Plot the OI vs wavelength.

    Parameters
    ----------
    dic : Dict
        Dictionary containing the data.
    key : str, optional
        The key in the dictionary to plot. The default is 'VIS2'.
    datatype : str, optional
        The datatype to plot. The default is 'VIS2'.
    showvis : bool, optional
        If True, the visibilities are plotted. The default is False.
    plot_errorbars : bool, optional
        If True, the errorbars are plotted. The default is True.
    correct_polynom : bool, optional
        If True, the polynomial is corrected. The default is False.
    timeVertOffset : int, optional
        The vertical offset of the time axis. The default is 0.
    stdevRange : bool, optional
        If True, the standard deviation is plotted. The default is False.
    normRange : bool, optional
        If True, the range is normalized. The default is False.
    stdevTime : bool, optional
        If True, the standard deviation is plotted. The default is False.
    """
    sta_index_cal = []
        
    # Get data from the input dictionary
    wl    = dic['WLEN'];
    nbWlen = len(wl)
    data  = dic[key][datatype];
    datae = dic[key][datatype + "ERR"];
    
    if correct_polynom:
        if normRange:
            l, h = normRange[0], normRange[1]
        else:
            l, h = 3, -3
        for i, j in enumerate(data):
            fit = np.polyfit(wl[l:h],data[i,l:h],correct_polynom)
            if key == 'VIS2' or key == 'TF2' or key == 'FLUX' or datatype == 'CFLUX':
                data[i,:] = data[i,:] / np.polyval(fit,wl)
            else:
                data[i,:] = data[i,:] - np.polyval(fit,wl)
    
    if key == 'FLUX':
        sta_index = dic[key]['STA_INDEX']
        sta_index_cal.append([sta_index])
    else:
        sta_index = np.sort(dic[key]['STA_INDEX'])
        sta_index_cal.append(sta_index)
                            
    # Get the unique station indices from the data
    sta_indices = np.unique(sta_index_cal, axis=0)
    if key == 'VIS' or key == 'VIS2' or key == 'TF2':
        n_max_config = np.nanmax([6, sta_indices.shape[0]])
        n_plot_rows = 3
    elif key == 'FLUX':
        n_max_config = 1#np.nanmax([4, sta_indices.shape[0]])
        n_plot_rows = 2
    elif key == 'T3':
        n_max_config = np.nanmax([4, sta_indices.shape[0]])
        n_plot_rows = 2
            
    print(n_max_config)
        
    fig1, axs1 = plt.subplots(n_plot_rows, 2, figsize=(15, 16), sharex=True, sharey=True)
    axs1 = axs1.ravel()
        
    #print datae.shape
    for i in range(int(len(data)//n_max_config)):
        label = datatype
        
        for j in range(n_max_config):
            idx = int(j*n_max_config+i);
            
            # Take square root of observable
            if showvis == True: 
                data[idx, :] = np.sqrt(data[idx, :])
                datae[idx, :] = 0.5 * datae[idx, :] / np.sqrt(data[idx, :])
                if key == 'VIS2' or key == 'TF2':
                    if key == 'VIS2':
                        label = 'VIS'
                    elif key == 'TF2':
                        label = 'TF'
                
            axs1[i].plot(wl * 1e6, data[idx, :]+timeVertOffset*j)
            axs1[i].set_ylabel(label)
            axs1[i].set_xlabel(r"$\lambda\ (\mu\mathrm{m}$)")
            if stdevRange:
                if datatype == 'CFLUX':
                    off = 1.05;
                else:
                    off = 0.5;
                axs1[i].text(wl[stdevRange[1]] * 1e6, timeVertOffset*j+off, "st. dev. in continuum"+str(round(np.std(data[idx, stdevRange[0]:stdevRange[1]]),3)))
                
            if plot_errorbars == True:
                axs1[i].errorbar(wl * 1e6, data[idx, :],
                    yerr=datae[idx, :],
                    ecolor = 'grey', alpha = 0.25, capsize = 0.5,elinewidth = 1)
            
        
        if stdevTime == True:
            dat = np.reshape(data,[int(len(data)//n_max_config),n_max_config,nbWlen])
            std = np.std(dat,2)
            if stdevRange:
                axs1[i].text(wl[stdevRange[1]] * 1e6, timeVertOffset*j+off, "st. dev. vs. time"+str(round(np.mean(np.std(dat[:,i, stdevRange[0]:stdevRange[1]],axis=0),3))))
                         
    if datatype == 'VIS2' or datatype == 'TF2' or datatype == 'CFLUX':
        plt.ylim([-0.1,1.1+timeVertOffset*(n_max_config-1)*1.1])
    else:
        try:
            mn = np.nanmin(data);
            if mn > 0:
                mn = mn*0.9;
            else:
                mn = mn*1.1
            plt.ylim([mn,np.nanmax(data+timeVertOffset*(n_max_config-1))*1.1])
            #plt.ylim([-30,60])
        except:
            pass
    if showvis == True:
        if datatype == 'VIS2':
            datatype = 'VIS'
        elif datatype == 'TF2':
            datatype = 'TF'
    fig1.suptitle(f'{datatype} vs. wavelength')
    plt.show()


###############################################################################
def show_oi_vs_time(list_of_dicts: Dict, wlenRange: List[float],
                    key: str = "VIS2", datatype: str = "VIS2",
                    showvis: Optional[bool] = False,
                    plot_errorbars: Optional[bool] = True) -> None:
    """This function shows the selected oifits data (flux, visibility, closure phase etc.) as a function of time.
    It reads data from multiple oifits files in a given directory.
    The data to be plotted can be filtered with the filter_oi_list function.

    Example usage:
    filtered_list_of_dicts = filter_oi_list(list_of_dicts,dates=["2018-03-14"],bands=['LM'],spectral_resolutions=['MED'],DIT_range=[0,0.2],targets=['l pup'])
    show_oi_vs_time(filtered_list_of_dicts, [3.5, 3.95], key="VIS2", datatype='VIS2') #[3.5, 3.95] [10.2,10.9]
    """
    # DicoFlux_L    
    # starFlux={'lamPyx': 28.105587, 'HD56618': 196.27484, 'omeLup': 115.3, 'Khya': 60.099998, 'AlfCMa': 995.79999, 'HD92397': 129.69469, 'epsCru': 217.05141, 'HD117597': 8.8213739, 'HD151843': 2.9975743, 'etaCen': 20.529898, 'delCen': 30.857565, 'HD100546': 20.5, 'HD107625': 0.49486896, 'BCen': 72.300003, 'HD162415': 12.079382, 'HD114873': 20.785559, 'HD103515': 24.602552, 'HD104329': 2.9723437, 'HD107869': 29.714039, 'muCen': 9.1099892, 'cVel': 108.7, 'tetCMa': 164.47496, 'HD109960': 17.799999, 'HD148984': 14.406607, 'HD150071': 0.31112576, 'HD105381': 12.466991, 'delSco': 32.341614, 'V1068Sco': 317.78168, 'HSco': 191.98087, 'HD140143': 9.6999998, 'HD97705': 0.043862239, 'HD148441': 0.99144912, '27Sco': 97.5, 'TetCen': 377.60001, 'HD102766': 145.40553, 'HD103975': 2.9723437, 'HD108110': 13.865284, 'HD97765': 0.12568563, 'HD113601': 5.0366859, 'HD148730': 0.20337902, '3C273': 0.059742939,'cVel': 108.7, 'omeLup': 115.3, 'Khya': 60.099998, 'AlfCMa': 995.79999, 'HD92397': 129.69469, 'epsCru': 217.05141, 'delCen': 30.857565, 'HSco': 191.98087, 'HD94495': 30.6, 'HD150071': 0.31112576, 'HD162189': 75.758179, 'BCen': 72.300003, 'TetCen': 377.60001, 'HD100546': 20.5, 'HD103515': 24.602552, 'HD97765': 0.12568563,'V532Car': 92.1,'epsSco': 394.2,'VV1068Sco': 100.,'upsLib': 200.,'gamAql': 600.}

    # DicoFluxN
    starFlux = {'lamPyx': 4.1238999, 'HD56618': 30.814779, 'omeLup': 17.6, 'Khya': 8.5500002, 'AlfCMa': 118.0, 'HD92397': 20.248549, 'epsCru': 32.011219, 'HD117597': 1.3875213, 'HD151843': 0.51499999, 'etaCen': 5.8867059, 'delCen': 14.54, 'HD100546': 45.694008, 'HD107625': 0.054777335, 'BCen': 10.616354, 'HD162415': 1.9299999, 'HD114873': 3.0, 'HD103515': 5.0599999, 'HD104329': 0.42240164, 'HD107869': 5.4688506, 'muCen': 2.6417792, 'cVel': 15.4, 'tetCMa': 30.499096, 'HD109960': 2.9200001, 'HD148984': 2.2, 'HD150071': 0.032316711, 'HD105381': 1.6773663, 'delSco': 12.018785, 'V1068Sco': 50.355038, 'HSco': 30.69574, 'HD140143': 1.71, 'HD97705': 0.0048251506, 'HD148441': 0.10740463, '27Sco': 10.9, 'TetCen': 48.400002, 'HD102766': 37.23, 'HD103975': 0.42240164, 'HD108110': 2.3499999, 'HD97765': 0.013920446, 'HD113601': 0.76300001, 'HD148730': 0.022324009, '3C273': 0.0067190719,'cVel': 15.4, 'omeLup': 17.6, 'Khya': 8.5500002, 'AlfCMa': 118.0, 'HD92397': 20.248549, 'epsCru': 32.011219, 'delCen': 14.54, 'HSco': 30.69574, 'HD94495': 5.02, 'HD150071': 0.032316711, 'HD162189': 13.1, 'BCen': 10.616354, 'TetCen': 48.400002, 'HD100546': 45.694008, 'HD103515': 5.0599999, 'HD97765': 0.013920446}
    #starFlux={'lamPyx': 16.899242, 'HD56618': 123.59395, 'omeLup': 64.288948, 'Khya': 34.413467, 'AlfCMa': 535.90002, 'HD92397': 79.279999, 'epsCru': 143.89999, 'HD117597': 5.3375411, 'HD151843': 1.8043216, 'etaCen': 18.671946, 'delCen': 27.200001, 'HD100546': 13.9, 'HD107625': 0.2566978, 'BCen': 40.465576, 'HD162415': 10.16, 'HD114873': 13.8, 'HD103515': 23.01, 'HD104329': 1.7560067, 'HD107869': 17.4, 'muCen': 6.5473495, 'cVel': 67.809998, 'tetCMa': 92.122192, 'HD109960': 9.8630896, 'HD148984': 5.6999998, 'HD150071': 0.16719383, 'HD105381': 5.1822686, 'delSco': 32.172329, 'V1068Sco': 198.0, 'HSco': 101.3, 'HD140143': 3.5664213, 'HD97705': 0.024131674, 'HD148441': 0.55915225, '27Sco': 48.699944, 'TetCen': 184.3, 'HD102766': 52.769211, 'HD103975': 1.7560067, 'HD108110': 8.4795961, 'HD97765': 0.068838015, 'HD113601': 4.2104497, 'HD148730': 0.11200868, '3C273': 0.03288243,'cVel': 67.809998, 'omeLup': 64.288948, 'Khya': 34.413467, 'AlfCMa': 535.90002, 'HD92397': 79.279999, 'epsCru': 143.89999, 'delCen': 27.200001, 'HSco': 101.3, 'HD94495': 12.719484, 'HD150071': 0.16719383, 'HD162189': 37.186783, 'BCen': 40.465576, 'TetCen': 184.3, 'HD100546': 13.9, 'HD103515': 23.01, 'HD97765': 0.068838015}
    # ('L', {'omeLup': 115.3, 'HD56618': 196.27484, 'HD148730': 0.20337902, 'Khya': 60.099998, 'AlfCMa': 995.79999, 'HD148984': 14.406607, 'HSco': 191.98087, 'epsCru': 217.05141, 'HD148441': 0.99144912, 'delCen': 30.857565, 'HD107869': 29.714039, 'HD104329': 2.9723437, 'HD105381': 12.466991, 'HD107625': 0.49486896, 'BCen': 72.300003, 'TetCen': 377.60001, 'HD103515': 24.602552, 'HD108110': 13.865284, 'HD140143': 9.6999998})
    # ('N', {'omeLup': 17.6, 'HD56618': 30.814779, 'HD148730': 0.022324009, 'Khya': 8.5500002, 'AlfCMa': 118.0, 'HD148984': 2.2, 'HSco': 30.69574, 'epsCru': 32.011219, 'HD148441': 0.10740463, 'delCen': 14.54, 'HD107869': 5.4688506, 'HD104329': 0.42240164, 'HD105381': 1.6773663, 'HD107625': 0.054777335, 'BCen': 10.616354, 'TetCen': 48.400002, 'HD103515': 5.0599999, 'HD108110': 2.3499999, 'HD140143': 1.71})
    # ('M', {'omeLup': 64.288948, 'HD56618': 123.59395, 'HD148730': 0.11200868, 'Khya': 34.413467, 'AlfCMa': 535.90002, 'HD148984': 5.6999998, 'HSco': 101.3, 'epsCru': 143.89999, 'HD148441': 0.55915225, 'delCen': 27.200001, 'HD107869': 17.4, 'HD104329': 1.7560067, 'HD105381': 5.1822686, 'HD107625': 0.2566978, 'BCen': 40.465576, 'TetCen': 184.3, 'HD103515': 23.01, 'HD108110': 8.4795961, 'HD140143': 3.5664213})

    # check if list is not empty:
    if list_of_dicts:
        target_names_cal = []
        MJD_arr_cal = []
        arr_cal = []
        err_arr_cal = []
        sta_index_cal = []

        target_names_sci = []
        MJD_arr_sci = []
        arr_sci = []
        err_arr_sci = []
        sta_index_sci = []

        for dic in list_of_dicts:
            wl = np.array(dic['WLEN'])
            wlenRange_idx = np.logical_and(wl > wlenRange[0] / 1.0e6, wl < wlenRange[1] / 1.0e6)
            category = dic['CATEGORY'].lower()
            if 'cal' in category:
                try:
                    datay    = np.array(dic[key][datatype])
                    datayerr = np.array(dic[key][datatype+'ERR'])
                    datax    = np.array(dic[key]["TIME"])
                    n_rows   = datay.shape[0]
                    # print datay.shape
                    for i in range(n_rows):
                        arr_cal.append(robust.mean(datay[i, wlenRange_idx]))
                        err_arr_cal.append(robust.mean(datayerr[i, wlenRange_idx]))
                        MJD_arr_cal.append(datax[i])
                        target_names_cal.append(dic['TARGET'])
                        if key == 'FLUX':
                            sta_index = dic[key]['STA_INDEX'][i]
                            sta_index_cal.append([sta_index])
                        else:
                            sta_index = np.sort(dic[key]['STA_INDEX'][i])
                            sta_index_cal.append(sta_index)
                except Exception:
                    print (dic['TARGET'], dic['DATEOBS'], "No CAL data found.")
            elif 'sci' in category:
                try:
                    datay    = np.array(dic[key][datatype])
                    datayerr = np.array(dic[key][datatype+'ERR'])
                    datax    = np.array(dic[key]["TIME"])
                    n_rows   = datay.shape[0]
                    # print datay.shape
                    for i in range(n_rows):
                        arr_sci.append(robust.mean(datay[i, wlenRange_idx]))
                        err_arr_sci.append(robust.mean(datayerr[i, wlenRange_idx]))
                        MJD_arr_sci.append(datax[i])
                        target_names_sci.append(dic['TARGET'])
                        if key == 'FLUX':
                            sta_index = dic[key]['STA_INDEX'][i]
                            sta_index_sci.append([sta_index])
                        else:
                            sta_index = np.sort(dic[key]['STA_INDEX'][i])
                            sta_index_sci.append(sta_index)
                except Exception:
                    print (dic['TARGET'], dic['DATEOBS'], "No SCI data found.")
            sta_names = dic['STA_NAME']

        target_names_cal = np.array(target_names_cal)
        MJD_arr_cal = np.array(MJD_arr_cal)
        arr_cal = np.array(arr_cal)
        err_arr_cal = np.array(err_arr_cal)
        sta_index_cal = np.array(sta_index_cal)

        target_names_sci = np.array(target_names_sci)
        MJD_arr_sci = np.array(MJD_arr_sci)
        arr_sci = np.array(arr_sci)
        err_arr_sci = np.array(err_arr_sci)
        sta_index_sci = np.array(sta_index_sci)

        sta_indices = np.unique(sta_index_cal, axis=0)
        if key == 'VIS' or key == 'VIS2' or key == 'TF2':
            n_max_config = np.nanmax([6, sta_indices.shape[0]])
            n_plot_rows = 3
        elif key == 'T3' or key == 'FLUX':
            n_max_config = np.nanmax([4, sta_indices.shape[0]])
            n_plot_rows = 2

        if len(MJD_arr_cal) > 0 and len(MJD_arr_sci) > 0:
            MJD_range = [np.nanmin([np.nanmin(MJD_arr_cal), np.nanmin(MJD_arr_sci)]),
                         np.nanmax([np.nanmax(MJD_arr_cal), np.nanmax(MJD_arr_sci)])]
        elif len(MJD_arr_sci) > 0:
            MJD_range = [np.nanmin(MJD_arr_sci), np.nanmax(MJD_arr_sci)]
        elif len(MJD_arr_cal) > 0:
            MJD_range = [np.nanmin(MJD_arr_cal), np.nanmax(MJD_arr_cal)]
        else:
            MJD_range = [0.0, 1.0]
        text_width_MJD = (MJD_range[1] - MJD_range[0]) / 20.0

        _, axs1 = plt.subplots(n_plot_rows, 2, figsize=(15, 16), sharex=True, sharey=True)
        axs1 = axs1.ravel()
        for i in range(n_max_config):
            if datatype == 'DPHI' or datatype == 'CLOS':
                axs1[i + 0].plot(MJD_range, [0.0, 0.0], '-', color='gray', lw=1.5)
            if len(sta_index_cal) > 0:
                idxst = np.all(sta_index_cal == sta_indices[i], axis=1)
                if len(arr_cal[idxst]) > 0:
                    label = datatype +' cal'
                    if plot_errorbars == True:
                        if showvis == True:
                            if key == 'VIS2' or key == 'TF2':
                                if key == 'VIS2':
                                    label = 'VIS' + ' cal'
                                elif key == 'TF2':
                                    label = 'TF' + ' cal'
                                axs1[i].errorbar(MJD_arr_cal[idxst], np.sqrt(arr_cal[idxst]),
                                                 yerr=0.5 * err_arr_cal[idxst] / np.sqrt(arr_cal[idxst]),
                                                 fmt='o', color='blue', elinewidth=1.0,
                                                 label=label)
                        else:
                            axs1[i].errorbar(MJD_arr_cal[idxst], arr_cal[idxst],
                                             yerr=err_arr_cal[idxst],
                                             fmt='o', color='blue', elinewidth=1.0,
                                             label=label)
                    else:
                        if showvis == True:
                            if key == 'VIS2' or key == 'TF2':
                                if key == 'VIS2':
                                    label = 'VIS' + ' sci'
                                elif key == 'TF2':
                                    label = 'TF' + ' sci'
                                axs1[i].errorbar(MJD_arr_cal[idxst], np.sqrt(arr_cal[idxst]),
                                                 fmt='o', color='blue', elinewidth=1.0,
                                                 label=label)
                        else:
                            axs1[i].errorbar(MJD_arr_cal[idxst], arr_cal[idxst],
                                             fmt='o', color='blue', elinewidth=1.0,
                                             label=label)
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if MJD_arr_cal[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_cal[idxst][j]):
                                ymin, ymax = axs1[i + 0].get_ylim()
                                axs1[i].text(MJD_arr_cal[idxst][j], ymax * 1.05,
                                             target_names_cal[idxst][j].replace('_', ' '), rotation=90,
                                             va='bottom')
                                text_tag_flag = 0
                                prev_text_MJD = MJD_arr_cal[idxst][j]
                                prev_target_name = target_names_cal[idxst][j]
            if len(sta_index_sci) > 0:
                label = datatype +' sci'
                idxst = np.all(sta_index_sci == sta_indices[i], axis=1)
                if len(arr_sci[idxst]) > 0:
                    if plot_errorbars == True:
                        if showvis == True:
                            if key == 'VIS2' or key == 'TF2':
                                if key == 'VIS2':
                                    label = 'VIS' + ' sci'
                                elif key == 'TF2':
                                    label = 'TF' + ' sci'
                                axs1[i].errorbar(MJD_arr_sci[idxst], np.sqrt(arr_sci[idxst]),
                                                 yerr=0.5 * err_arr_sci[idxst] / np.sqrt(arr_sci[idxst]),
                                                 fmt='o', color='red', elinewidth=1.0,
                                                 label=label)
                        else:
                            axs1[i].errorbar(MJD_arr_sci[idxst], arr_sci[idxst], yerr=err_arr_sci[idxst],
                                                 fmt='o', color='red', elinewidth=1.0,
                                                 label=label)
                    else:
                        if showvis == True:
                            if key == 'VIS2' or key == 'TF2':
                                if key == 'VIS2':
                                    label = 'VIS' + ' sci'
                                elif key == 'TF2':
                                    label = 'TF' + ' sci'
                                axs1[i].errorbar(MJD_arr_sci[idxst], np.sqrt(arr_sci[idxst]),
                                                 fmt='o', color='red', elinewidth=1.0,
                                                 label=label)
                        else:
                            axs1[i].errorbar(MJD_arr_sci[idxst], arr_sci[idxst],
                                         fmt='o', color='red', elinewidth=1.0,
                                         label=label)
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if MJD_arr_sci[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_sci[idxst][j]):
                                ymin, ymax = axs1[i + 0].get_ylim()
                                axs1[i].text(MJD_arr_sci[idxst][j], ymax*1.05, target_names_sci[idxst][j].replace('_', ' '),
                                             rotation=90, va='bottom', color='darkred')
                                text_tag_flag = 0
                                prev_text_MJD = MJD_arr_sci[idxst][j]
                                prev_target_name = target_names_sci[idxst][j]
            if key == 'VIS' or key == 'VIS2' or key == 'TF2':
                axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                      sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]
            elif key == 'T3':
                axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                      sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0] + ' - ' + \
                      sta_names[sta_indices[i, 2] == dic['STA_INDEX']][0]
            elif key == 'FLUX':
                axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0]
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            axs1[i].text(0.05, 0.95, axlabel, horizontalalignment='left', verticalalignment='top',
                         transform=axs1[i].transAxes, bbox=props)
            if i == 0:
                leg = axs1[i].legend(loc='upper right')
                leg.get_frame().set_alpha(0.5)
            if datatype == 'VIS2' or datatype == 'TF2':
                axs1[i].set_ylim([-0.1, 1.1])
            else:
                try:
                    axs1[i].set_ylim([np.nanmin([np.nanmin(arr_cal),np.nanmin(arr_sci)]), np.nanmax([np.nanmax(arr_cal),np.nanmax(arr_sci)])])
                except:
                    pass
            ylabel = datatype
            if showvis:
                if datatype == 'VIS2':
                    ylabel = 'VIS'
                    datatype = 'VIS'
                elif datatype == 'TF2':
                    ylabel = 'TF'
                    datatype = 'TF'
            axs1[i].set_ylabel(ylabel)
            axs1[i].set_xlabel(r'$\mathrm{MJD}$')
        plt.suptitle(r'$\mathrm{' + datatype + r'\ vs.\ time}$')
        plt.tight_layout()
        plt.show()


def show_vis2_tf2_vs_time(
        list_of_dicts: List[Dict], wlenRange: List[float],
        showvis: Optional[bool] = False,
        saveplots: Optional[bool] = False,
        output_path: Optional[Path] = None,
        plot_errorbars: Optional[bool] = True,
        magic_numbers: Optional[bool] = False, **kwargs) -> None:
    """Plot the vis2 and tf2 vs time.

    Parameters
    ----------
    list_of_dicts : list of dict
        List of oifits-files.
    wlenRange : list
        Wavelength range for the plot.
    showvis : bool, optional
        If True, plot visibilities (V) instead of V^2: V is calculated from V^2 (not from the VISAMP table). The default is False.
    saveplots : bool, optional
        If True, the plots are saved in the output_path. The default is False.
    output_path : pathlib.Path, optional
        Path to save the plots. The default is None.
    plot_errorbars : bool, optional
        If True, plot errorbars. The default is True.
    magic_numbers : bool, optional
        If True, plot the magic numbers. The default is False.
    """
    print('Go !')
    #starFlux={'Teta Pyx':27, 'psi Vir':44,'tet Cen':48, \
    #      'del Vir':160,'eps Sco':52,'del Oph':110, \
    #      'FZ Lib':17,'del Sco':15,'HD171094':34, \
    #      'e Aql':27,'29 Cap':31,'alf Ant':19, \
    #      'HD99333':20,'C01 Cen':22,'alpha Arae':10, \
    #      'HD138505':18,'HD142527':12,'HD161560':0.4, \
    #      'HD171094':34,'AV Mic':18,'HD181041':6,  \
    #      'HD138492':8,'HD147929':8,'HD148255':3, \
    #      'G Sco':26,'HD156936':3,'HD161849':8, \
    #      'HD165413':5,'HD186765':9,'Mu Hya':30, \
    #      'CH Vir':20,'HD 126111':6,'LY TrA':7, \
    #      'ET Vir':29,'H Sco':30,'HD177507':17, \
    #      'V345 Tel':9,'RW Lup':18,'HD 150798':144, \
    #      'BM Sco':85,'RX Tel':82,'HD328913':100, \
    #      'Eps Oph':16,'HD182669':2.5,'Nu Hya':30}    # check if list is not empty:

    #DicoFlux_L    
    #starFlux={'lamPyx': 28.105587, 'HD56618': 196.27484, 'omeLup': 115.3, 'Khya': 60.099998, 'AlfCMa': 995.79999, 'HD92397': 129.69469, 'epsCru': 217.05141, 'HD117597': 8.8213739, 'HD151843': 2.9975743, 'etaCen': 20.529898, 'delCen': 30.857565, 'HD100546': 20.5, 'HD107625': 0.49486896, 'BCen': 72.300003, 'HD162415': 12.079382, 'HD114873': 20.785559, 'HD103515': 24.602552, 'HD104329': 2.9723437, 'HD107869': 29.714039, 'muCen': 9.1099892, 'cVel': 108.7, 'tetCMa': 164.47496, 'HD109960': 17.799999, 'HD148984': 14.406607, 'HD150071': 0.31112576, 'HD105381': 12.466991, 'delSco': 32.341614, 'V1068Sco': 317.78168, 'HSco': 191.98087, 'HD140143': 9.6999998, 'HD97705': 0.043862239, 'HD148441': 0.99144912, '27Sco': 97.5, 'TetCen': 377.60001, 'HD102766': 145.40553, 'HD103975': 2.9723437, 'HD108110': 13.865284, 'HD97765': 0.12568563, 'HD113601': 5.0366859, 'HD148730': 0.20337902, '3C273': 0.059742939,'cVel': 108.7, 'omeLup': 115.3, 'Khya': 60.099998, 'AlfCMa': 995.79999, 'HD92397': 129.69469, 'epsCru': 217.05141, 'delCen': 30.857565, 'HSco': 191.98087, 'HD94495': 30.6, 'HD150071': 0.31112576, 'HD162189': 75.758179, 'BCen': 72.300003, 'TetCen': 377.60001, 'HD100546': 20.5, 'HD103515': 24.602552, 'HD97765': 0.12568563,'V532Car': 92.1,'epsSco': 394.2,'VV1068Sco': 100.,'upsLib': 200.,'gamAql': 600.}

    # dicoFluxL={'lamPyx': 28.105587, 'HD56618': 196.27484, 'omeLup': 115.3, 'Khya': 60.099998, 'AlfCMa': 995.79999, 'HD92397': 129.69469, 'epsCru': 217.05141, 'HD117597': 8.8213739, 'HD151843': 2.9975743, 'etaCen': 20.529898, 'delCen': 30.857565, 'HD100546': 20.5, 'HD107625': 0.49486896, 'BCen': 72.300003, 'HD162415': 12.079382, 'HD114873': 20.785559, 'HD103515': 24.602552, 'HD104329': 2.9723437, 'HD107869': 29.714039, 'muCen': 9.1099892, 'cVel': 108.7, 'tetCMa': 164.47496, 'HD109960': 17.799999, 'HD148984': 14.406607, 'HD150071': 0.31112576, 'HD105381': 12.466991, 'delSco': 32.341614, 'V1068Sco': 317.78168, 'HSco': 191.98087, 'HD140143': 9.6999998, 'HD97705': 0.043862239, 'HD148441': 0.99144912, '27Sco': 97.5, 'TetCen': 377.60001, 'HD102766': 145.40553, 'HD103975': 2.9723437, 'HD108110': 13.865284, 'HD97765': 0.12568563, 'HD113601': 5.0366859, 'HD148730': 0.20337902, '3C273': 0.059742939,'cVel': 108.7, 'omeLup': 115.3, 'Khya': 60.099998, 'AlfCMa': 995.79999, 'HD92397': 129.69469, 'epsCru': 217.05141, 'delCen': 30.857565, 'HSco': 191.98087, 'HD94495': 30.6, 'HD150071': 0.31112576, 'HD162189': 75.758179, 'BCen': 72.300003, 'TetCen': 377.60001, 'HD100546': 20.5, 'HD103515': 24.602552, 'HD97765': 0.12568563}
    # dicoFluxN
    starFlux= {'lamPyx': 4.1238999, 'HD56618': 30.814779, 'omeLup': 17.6, 'Khya': 8.5500002, 'AlfCMa': 118.0, 'HD92397': 20.248549, 'epsCru': 32.011219, 'HD117597': 1.3875213, 'HD151843': 0.51499999, 'etaCen': 5.8867059, 'delCen': 14.54, 'HD100546': 45.694008, 'HD107625': 0.054777335, 'BCen': 10.616354, 'HD162415': 1.9299999, 'HD114873': 3.0, 'HD103515': 5.0599999, 'HD104329': 0.42240164, 'HD107869': 5.4688506, 'muCen': 2.6417792, 'cVel': 15.4, 'tetCMa': 30.499096, 'HD109960': 2.9200001, 'HD148984': 2.2, 'HD150071': 0.032316711, 'HD105381': 1.6773663, 'delSco': 12.018785, 'V1068Sco': 50.355038, 'HSco': 30.69574, 'HD140143': 1.71, 'HD97705': 0.0048251506, 'HD148441': 0.10740463, '27Sco': 10.9, 'TetCen': 48.400002, 'HD102766': 37.23, 'HD103975': 0.42240164, 'HD108110': 2.3499999, 'HD97765': 0.013920446, 'HD113601': 0.76300001, 'HD148730': 0.022324009, '3C273': 0.0067190719,'cVel': 15.4, 'omeLup': 17.6, 'Khya': 8.5500002, 'AlfCMa': 118.0, 'HD92397': 20.248549, 'epsCru': 32.011219, 'delCen': 14.54, 'HSco': 30.69574, 'HD94495': 5.02, 'HD150071': 0.032316711, 'HD162189': 13.1, 'BCen': 10.616354, 'TetCen': 48.400002, 'HD100546': 45.694008, 'HD103515': 5.0599999, 'HD97765': 0.013920446}
    # starFlux={'lamPyx': 16.899242, 'HD56618': 123.59395, 'omeLup': 64.288948, 'Khya': 34.413467, 'AlfCMa': 535.90002, 'HD92397': 79.279999, 'epsCru': 143.89999, 'HD117597': 5.3375411, 'HD151843': 1.8043216, 'etaCen': 18.671946, 'delCen': 27.200001, 'HD100546': 13.9, 'HD107625': 0.2566978, 'BCen': 40.465576, 'HD162415': 10.16, 'HD114873': 13.8, 'HD103515': 23.01, 'HD104329': 1.7560067, 'HD107869': 17.4, 'muCen': 6.5473495, 'cVel': 67.809998, 'tetCMa': 92.122192, 'HD109960': 9.8630896, 'HD148984': 5.6999998, 'HD150071': 0.16719383, 'HD105381': 5.1822686, 'delSco': 32.172329, 'V1068Sco': 198.0, 'HSco': 101.3, 'HD140143': 3.5664213, 'HD97705': 0.024131674, 'HD148441': 0.55915225, '27Sco': 48.699944, 'TetCen': 184.3, 'HD102766': 52.769211, 'HD103975': 1.7560067, 'HD108110': 8.4795961, 'HD97765': 0.068838015, 'HD113601': 4.2104497, 'HD148730': 0.11200868, '3C273': 0.03288243,'cVel': 67.809998, 'omeLup': 64.288948, 'Khya': 34.413467, 'AlfCMa': 535.90002, 'HD92397': 79.279999, 'epsCru': 143.89999, 'delCen': 27.200001, 'HSco': 101.3, 'HD94495': 12.719484, 'HD150071': 0.16719383, 'HD162189': 37.186783, 'BCen': 40.465576, 'TetCen': 184.3, 'HD100546': 13.9, 'HD103515': 23.01, 'HD97765': 0.068838015}

    # ('L', {'omeLup': 115.3, 'HD56618': 196.27484, 'HD148730': 0.20337902, 'Khya': 60.099998, 'AlfCMa': 995.79999, 'HD148984': 14.406607, 'HSco': 191.98087, 'epsCru': 217.05141, 'HD148441': 0.99144912, 'delCen': 30.857565, 'HD107869': 29.714039, 'HD104329': 2.9723437, 'HD105381': 12.466991, 'HD107625': 0.49486896, 'BCen': 72.300003, 'TetCen': 377.60001, 'HD103515': 24.602552, 'HD108110': 13.865284, 'HD140143': 9.6999998})
    # ('N', {'omeLup': 17.6, 'HD56618': 30.814779, 'HD148730': 0.022324009, 'Khya': 8.5500002, 'AlfCMa': 118.0, 'HD148984': 2.2, 'HSco': 30.69574, 'epsCru': 32.011219, 'HD148441': 0.10740463, 'delCen': 14.54, 'HD107869': 5.4688506, 'HD104329': 0.42240164, 'HD105381': 1.6773663, 'HD107625': 0.054777335, 'BCen': 10.616354, 'TetCen': 48.400002, 'HD103515': 5.0599999, 'HD108110': 2.3499999, 'HD140143': 1.71})
    # ('M', {'omeLup': 64.288948, 'HD56618': 123.59395, 'HD148730': 0.11200868, 'Khya': 34.413467, 'AlfCMa': 535.90002, 'HD148984': 5.6999998, 'HSco': 101.3, 'epsCru': 143.89999, 'HD148441': 0.55915225, 'delCen': 27.200001, 'HD107869': 17.4, 'HD104329': 1.7560067, 'HD105381': 5.1822686, 'HD107625': 0.2566978, 'BCen': 40.465576, 'TetCen': 184.3, 'HD103515': 23.01, 'HD108110': 8.4795961, 'HD140143': 3.5664213})

    if list_of_dicts:
        # colors: BCD: out-out, in-in, in-out, out-in
        BCD_configs = np.array([[0, 0], [1, 1], [1, 0], [0, 1]])
        BCD_labels = np.array(['OUT-OUT', 'IN-IN', 'IN-OUT', 'OUT-IN'])
        BCD_markers = np.array(['o', 's', 'd', 'p'])
        V2_cal_colors = np.array(['darkseagreen', 'yellowgreen', 'olivedrab', 'darkkhaki'])
        V2_colors = np.array(['red', 'orange', 'salmon', 'pink'])
        TF2_colors = np.array(['blue', 'dodgerblue', 'royalblue', 'indigo'])

        target_names_cal = []
        target_fluxes_cal = []
        V2_BCD_arr_cal = []
        V2_MJD_arr_cal = []
        V2_arr_cal = []
        V2err_arr_cal = []
        V2_sta_index_cal = []

        target_names_CP_cal = []
        target_fluxes_CP_cal = []
        CP_BCD_arr_cal = []
        CP_MJD_arr_cal = []
        CP_arr_cal = []
        CPerr_arr_cal = []
        CP_sta_index_cal = []

        target_names_TF2 = []
        target_fluxes_TF2 = []
        TF2_BCD_arr = []
        TF2_MJD_arr = []
        TF2_arr = []
        TF2err_arr = []
        TF_arr = []
        TFerr_arr = []
        TF2_sta_index = []

        target_names = []
        target_fluxes = []
        V2_BCD_arr = []
        V2_MJD_arr = []
        V2_arr = []
        V2err_arr = []
        V2_sta_index = []

        target_names_CP = []
        target_fluxes_CP = []
        CP_BCD_arr = []
        CP_MJD_arr = []
        CP_arr = []
        CPerr_arr = []
        CP_sta_index = []

        #print(list_of_dicts[0])
        for dic in list_of_dicts:
            print(dic['TARGET'])
            wl = np.array(dic['WLEN'])
            nwave=np.shape(wl)[0]
            # magic_a_L=np.array([17855.33103,1585.766286,-442130.2517,-84923.15308,32355,58401,178739.575])
            # magic_b_L=np.array([0.933883878,1.001444557,3.112701998,1.440708905,0.797069729,0.05296366])
            # magic_a_M=np.array([40981.87053,-12155.63552,-17532.4504,80916.31496,22219.88766,-84838.5665])
            # magic_b_M=np.array([0.809669717,1.081702569,1.330259383,0.648748387,0.883684074,1.225873326])
            magic_a_L=np.array([17855.33103,1585.766286,178739.575,32355.58401,-84923.15308,-442130.2517])
            magic_b_L=np.array([0.933883878,1.001444557,0.05296366,0.797069729,1.440708905,3.112701998])
            magic_a_M=np.array([40981.87053,-12155.63552,-84838.5665,22219.88766,80916.31496,-17532.4504])
            magic_b_M=np.array([0.809669717,1.081702569,1.225873326,0.883684074,0.648748387,1.330259383])
            ind_L=np.where(wl < 4.2e-6)
            ind_M=np.where(wl > 4.2e-6)
            print('wl = ',wl[ind_L])
            magic_curve_V2=np.zeros((6,nwave),dtype=np.float64)
            for i in range(6):
                magic_curve_V2[i,ind_L]=wl[ind_L]*magic_a_L[i]+magic_b_L[i]
                magic_curve_V2[i,ind_M]=wl[ind_M]*magic_a_M[i]+magic_b_M[i]
            wlenRange_idx = np.logical_and(wl > wlenRange[0] / 1.0e6, wl < wlenRange[1] / 1.0e6)
            if sum(wlenRange_idx) > 0:
                category = dic['CATEGORY'].lower()
                print(category)
                if 'cal' in category:
                    try:
                        datay = np.array(dic['VIS2']['VIS2'])
                        datayerr = np.array(dic['VIS2']['VIS2ERR'])
                        datax = np.array(dic['VIS2']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):   
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            V2_arr_cal.append(robust.mean(datay[i, wlenRange_idx]))
                            V2err_arr_cal.append(robust.mean(datayerr[i, wlenRange_idx]))
                            V2_BCD_arr_cal.append([BCD1, BCD2])
                            V2_MJD_arr_cal.append(datax[i])
                            target_names_cal.append(dic['TARGET'])
                            target_fluxes_cal.append(dic['TARGET_FLUX'])
                            sta_index = np.sort(dic['VIS2']['STA_INDEX'][i])
                            V2_sta_index_cal.append(sta_index)
                    except Exception:
                        print (dic['TARGET'], dic['DATEOBS'], "No CAL VIS2 data found.")
                    try:
                        datay = np.array(dic['T3']['CLOS'])
                        datayerr = np.array(dic['T3']['CLOSERR'])
                        datax = np.array(dic['T3']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            CP_arr_cal.append(robust.mean(datay[i, wlenRange_idx]))
                            CPerr_arr_cal.append(robust.mean(datayerr[i, wlenRange_idx]))
                            CP_BCD_arr_cal.append([BCD1, BCD2])
                            CP_MJD_arr_cal.append(datax[i])
                            target_names_CP_cal.append(dic['TARGET'])
                            target_fluxes_CP_cal.append(dic['TARGET_FLUX'])
                            sta_index = np.sort(dic['T3']['STA_INDEX'][i])
                            CP_sta_index_cal.append(sta_index)
                    except Exception:
                        print (dic['TARGET'], dic['DATEOBS'], "No CAL CP data found.")
                    try:
                        print(np.shape(dic['TF2']['STA_INDEX']))
                        datay = np.array(dic['TF2']['TF2'])
                        datayerr = np.array(dic['TF2']['TF2ERR'])
                        datax = np.array(dic['TF2']["TIME"])
                        n_rows = datay.shape[0]
                        nB=6
                        nexp=np.int(n_rows/6)
                        
                        if (magic_numbers == True and dic['BCD1NAME'] == 'IN' and dic['BCD2NAME'] == 'IN'):
                            print('sta_index = ',dic['TF2']['STA_INDEX'])
                            for i in range(nexp):
                                k=i*nB
                                datay[k:k+nB,:]=datay[k:k+nB,:]*magic_curve_V2[0:nB,:]

                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            TF2_arr.append(robust.mean(datay[i, wlenRange_idx]))
                            TF2err_arr.append(np.mean(datayerr[i, wlenRange_idx]))
                            TF2_BCD_arr.append([BCD1, BCD2])
                            TF2_MJD_arr.append(datax[i])
                            target_names_TF2.append(dic['TARGET'])
                            target_fluxes_TF2.append(dic['TARGET_FLUX'])
                            sta_index = np.sort(dic['TF2']['STA_INDEX'][i])
                            print('sta_index sorted = ',np.sort(dic['TF2']['STA_INDEX'][i]))
                            TF2_sta_index.append(sta_index)
                    except NameError:
                        print(dic['TARGET'], dic['DATEOBS'], "No CAL TF2 data found.")
                if 'sci' in category:
                    try:
                        datay = np.array(dic['VIS2']['VIS2'])
                        datayerr = np.array(dic['VIS2']['VIS2ERR'])
                        datax = np.array(dic['VIS2']["TIME"])
                        n_rows = datay.shape[0]
                        nB=6
                        nexp=np.int(n_rows/6)
                        if (magic_numbers == True and dic['BCD1NAME'] == 'IN' and dic['BCD2NAME'] == 'IN'):
                            for i in range(nexp):
                                k=i*nB
                                datay[k:k+nB,:]=datay[k:k+nB,:]*magic_curve_V2[0:nB,:]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            V2_arr.append(robust.mean(datay[i, wlenRange_idx]))
                            V2err_arr.append(robust.mean(datayerr[i, wlenRange_idx]))
                            V2_BCD_arr.append([BCD1, BCD2])
                            V2_MJD_arr.append(datax[i])
                            target_names.append(dic['TARGET'])
                            target_fluxes.append(dic['TARGET_FLUX'])
                            sta_index = np.sort(dic['VIS2']['STA_INDEX'][i])
                    except Exception:
                        print (dic['TARGET'], dic['DATEOBS'], "No SCI VIS2 data found.")
                    try:
                        datay = np.array(dic['T3']['CLOS'])
                        datayerr = np.array(dic['T3']['CLOSERR'])
                        datax = np.array(dic['T3']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            CP_arr.append(robust.mean(datay[i, wlenRange_idx]))
                            CPerr_arr.append(robust.mean(datayerr[i, wlenRange_idx]))
                            target_names_CP.append(dic['TARGET'])
                            target_fluxes_CP.append(dic['TARGET_FLUX'])
                            sta_index = np.sort(dic['T3']['STA_INDEX'][i])
                            CP_sta_index.append(sta_index)
                            CP_BCD_arr.append([BCD1, BCD2])
                            CP_MJD_arr.append(datax[i])
                    except Exception:
                        print (dic['TARGET'], dic['DATEOBS'], "No SCI CP data found.")
            else:
                print("Wavelength out of range.")

        sta_names = dic['STA_NAME']
        
        target_names_cal = np.array(target_names_cal)
        target_fluxes_cal = np.array(target_fluxes_cal)
        V2_BCD_arr_cal = np.array(V2_BCD_arr_cal)
        V2_MJD_arr_cal = np.array(V2_MJD_arr_cal)
        V2_arr_cal = np.array(V2_arr_cal)
        V2err_arr_cal = np.array(V2err_arr_cal)
        V2_sta_index_cal = np.array(V2_sta_index_cal)

        target_names_CP_cal = np.array(target_names_CP_cal)
        target_fluxes_CP_cal = np.array(target_fluxes_CP_cal)
        CP_BCD_arr_cal = np.array(CP_BCD_arr_cal)
        CP_MJD_arr_cal = np.array(CP_MJD_arr_cal)
        CP_arr_cal = np.array(CP_arr_cal)
        CPerr_arr_cal = np.array(CPerr_arr_cal)
        CP_sta_index_cal = np.array(CP_sta_index_cal)

        target_names_TF2 = np.array(target_names_TF2)
        target_fluxes_TF2 = np.array(target_fluxes_TF2)
        TF2_BCD_arr = np.array(TF2_BCD_arr)
        TF2_MJD_arr = np.array(TF2_MJD_arr)
        TF2_arr = np.array(TF2_arr)
        TF2err_arr = np.array(TF2err_arr)
        TF2_sta_index = np.array(TF2_sta_index)

        target_names = np.array(target_names)
        target_fluxes = np.array(target_fluxes)
        V2_BCD_arr = np.array(V2_BCD_arr)
        V2_MJD_arr = np.array(V2_MJD_arr)
        V2_arr = np.array(V2_arr)
        V2err_arr = np.array(V2err_arr)
        V2_sta_index = np.array(V2_sta_index)

        target_names_CP = np.array(target_names_CP)
        target_fluxes_CP = np.array(target_fluxes_CP)
        CP_BCD_arr = np.array(CP_BCD_arr)
        CP_MJD_arr = np.array(CP_MJD_arr)
        CP_arr = np.array(CP_arr)
        CPerr_arr = np.array(CPerr_arr)
        CP_sta_index = np.array(CP_sta_index)

        if len(V2_sta_index_cal) > 0:
            sta_indices = np.unique(V2_sta_index_cal, axis=0)
        elif len(V2_sta_index) > 0:
            sta_indices = np.unique(V2_sta_index, axis=0)
        else:
            print ("Data arrays empty. Quitting.")
            return
        n_max_config = np.nanmax([6, sta_indices.shape[0]])

        if len(V2_MJD_arr_cal) > 0 and len(V2_MJD_arr) > 0:
            MJD_range = [np.nanmin([np.nanmin(V2_MJD_arr_cal), np.nanmin(V2_MJD_arr)]),
                         np.nanmax([np.nanmax(V2_MJD_arr_cal), np.nanmax(V2_MJD_arr)])]
        elif len(V2_MJD_arr) > 0:
            MJD_range = [np.nanmin(V2_MJD_arr), np.nanmax(V2_MJD_arr)]
        elif len(V2_MJD_arr_cal) > 0:
            MJD_range = [np.nanmin(V2_MJD_arr_cal), np.nanmax(V2_MJD_arr_cal)]
        else:
            MJD_range = [0.0, 1.0]
        text_width_MJD = (MJD_range[1] - MJD_range[0]) / 20.0
        fig1, axs1 = plt.subplots(3, 2, figsize=(16, 16), sharex=True, sharey=True)
        axs1 = axs1.ravel()
        text_y = 1.02
        for i in range(n_max_config):
            if len(V2_sta_index_cal) > 0:
                idxst = np.all(V2_sta_index_cal == sta_indices[i], axis=1)
                if len(V2_arr_cal[idxst]) > 0:
                    if showvis:
                        label = 'V cal '
                    else:
                        label = 'V2 cal '
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(V2_BCD_arr_cal == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)

                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if V2_MJD_arr_cal[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_cal[idxst][j]):
                                axs1[i].text(V2_MJD_arr_cal[idxst][j], text_y, \
                                                 target_names_cal[idxst][j].replace('_', ' ')+ \
                                                 " ("+np.str(np.round(target_fluxes_cal[idxst][j],1))+"Jy)", rotation=90, \
                                                 va='bottom',fontsize=7)
                                text_tag_flag = 0
                                prev_text_MJD = V2_MJD_arr_cal[idxst][j]
                                prev_target_name = target_names_cal[idxst][j]

            if len(TF2_sta_index) > 0:
                if showvis:
                    label = 'TF '
                else:
                    label = 'TF2 '
                idxst = np.all(TF2_sta_index == sta_indices[i], axis=1)
                print('np.shape(idxst) = {0}'.format(np.shape(idxst)))
                print('np.shape(TF2_arr) = {0}'.format(np.shape(TF2_arr)))
                if len(TF2_arr[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(TF2_BCD_arr == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(TF2_arr[cidxst]) > 0:
                            if plot_errorbars:
                                if showvis:
                                    axs1[i].errorbar((TF2_MJD_arr[cidxst]-np.min(TF2_MJD_arr[cidxst]))*24., np.sqrt(np.abs(TF2_arr[cidxst])),
                                                     yerr=0.5 * TF2err_arr[cidxst] / np.sqrt(np.abs(TF2_arr[cidxst])),
                                                     fmt=BCD_markers[j], color=TF2_colors[j], elinewidth=3.5)

                                    print('Direct RMS={0},{1},{2}'.format(np.std(np.sqrt(np.abs(TF2_arr[cidxst]))),BCD_labels[j],sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]))
                                else:
                                    axs1[i].errorbar(TF2_MJD_arr[cidxst], TF2_arr[cidxst], yerr=TF2err_arr[cidxst],
                                                     fmt=BCD_markers[j], color=TF2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                            else:
                                if showvis:
                                    axs1[i].errorbar(TF2_MJD_arr[cidxst], np.sqrt(np.abs(TF2_arr[cidxst])),
                                                     fmt=BCD_markers[j], color=TF2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                else:
                                    axs1[i].errorbar(TF2_MJD_arr[cidxst], TF2_arr[cidxst],
                                                     fmt=BCD_markers[j], color=TF2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
            if len(V2_sta_index) > 0:
                if showvis:
                    label = 'V sci '
                else:
                    label = 'V2 sci '
                idxst = np.all(V2_sta_index == sta_indices[i], axis=1)
                if len(V2_arr[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(V2_BCD_arr == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(V2_arr[cidxst]) > 0:
                            if plot_errorbars:
                                if showvis:
                                    axs1[i].errorbar(V2_MJD_arr[cidxst], np.sqrt(V2_arr[cidxst]),
                                                     yerr=0.5 * V2err_arr[cidxst] / np.sqrt(V2_arr[cidxst]),
                                                     fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                    print('Direct RMS={0},{1},{2}'.format(np.std(np.sqrt(np.abs(V2_arr[cidxst]))),BCD_labels[j],sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]))
                                else:
                                    axs1[i].errorbar(V2_MJD_arr[cidxst], V2_arr[cidxst], yerr=V2err_arr[cidxst],
                                                     fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                    print('Direct RMS={0},{1},{2}'.format(np.std(np.abs(V2_arr[cidxst])),BCD_labels[j],sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]))
                            else:
                                if showvis:
                                    axs1[i].errorbar(V2_MJD_arr[cidxst], np.sqrt(V2_arr[cidxst]),
                                                     fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                else:
                                    axs1[i].errorbar(V2_MJD_arr[cidxst], V2_arr[cidxst],
                                                     fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if V2_MJD_arr[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names[idxst][j]):
                                axs1[i].text(V2_MJD_arr[idxst][j], text_y, target_names[idxst][j].replace('_', ' '),
                                             rotation=90, va='bottom', color='darkred',fontsize=9)
                                text_tag_flag = 0
                                prev_text_MJD = V2_MJD_arr[idxst][j]
                                prev_target_name = target_names[idxst][j]


            axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                      sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            axs1[i].text(0.05, 0.9, axlabel, horizontalalignment='left', verticalalignment='top',
                         transform=axs1[i].transAxes, bbox=props,fontsize=15)
            axs1[i].tick_params(axis='both',labelsize='16')
            axs1[i].set_ylim([0.0, 1.0])
            if showvis:
                ylabel = 'V'
            else:
                ylabel = '$V^2$'
            axs1[i].set_ylabel(ylabel,fontsize=17)
            axs1[i].set_xlabel('Time [MJD]',fontsize=17)
        if showvis:
            plt.suptitle(r'$V\mathrm{\ vs.\ time}$')
        else:
            plt.suptitle(r'$V^2\mathrm{\ vs.\ time}$')
        fig1.subplots_adjust(hspace=0.3, wspace=0.3)
        for i in range(4):
            plt.setp(axs1[i].get_xticklabels(), visible=True)
            x_axis = axs1[i].axes.get_xaxis()
            x_axis.get_label().set_visible(True)
        for i in range(1, 6):
            plt.setp(axs1[i].get_yticklabels(), visible=True)
            y_axis = axs1[i].axes.get_yaxis()
            y_axis.get_label().set_visible(True)
        if saveplots:
            if showvis:
                label = '_VIS_TF'
            else:
                label = '_VIS2_TF2'
            fig1.savefig(output_path / f"{label}.png", dpi=300)
            fig1.savefig(output_path / f"{label}.eps", format='eps', dpi=300)
            plt.close(fig1)

        if len(CP_sta_index_cal) > 0:
            CP_sta_indices = np.unique(CP_sta_index_cal, axis=0)
        else:
            CP_sta_indices = np.unique(CP_sta_index, axis=0)

        # print CP_sta_indices
        n_max_config = np.nanmax([4, CP_sta_indices.shape[0]])

        fig2, axs = plt.subplots(2, 2, figsize=(14, 12), sharex=True, sharey=True)
        axs = axs.ravel()
        text_y = 60

        for i in range(n_max_config):
            axs[i + 0].plot(MJD_range, [0.0, 0.0], '-', color='gray', lw=1.5)
            if len(CP_sta_index_cal) > 0:
                idxst = np.all(CP_sta_index_cal == CP_sta_indices[i], axis=1)
                print(np.shape(CP_MJD_arr_cal))
                print(np.shape(CP_arr_cal))
                if len(CP_arr_cal[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(CP_BCD_arr_cal == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(CP_arr_cal[cidxst]) > 0:
                            if plot_errorbars == True:
                                axs[i + 0].errorbar(CP_MJD_arr_cal[cidxst], CP_arr_cal[cidxst], yerr=CPerr_arr_cal[cidxst],
                                                fmt=BCD_markers[j], color=V2_cal_colors[j], elinewidth=1.5,
                                                label='CP cal ' + BCD_labels[j])
                                axs[i+0].set_ylim(-10,10)
                                print('Direct RMS={0},{1},{2}'.format(np.std(CP_arr_cal[cidxst]),BCD_labels[j],sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]))
                            else:
                                axs[i + 0].errorbar(CP_MJD_arr_cal[cidxst], CP_arr_cal[cidxst],
                                                    fmt=BCD_markers[j], color=V2_cal_colors[j], elinewidth=1.5,
                                                    label='CP cal ' + BCD_labels[j])
                                axs[i+0].set_ylim(-10,10)

                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if CP_MJD_arr_cal[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_CP_cal[idxst][j]):
                                ymin, ymax = axs[i + 0].get_ylim()
                                axs[i + 0].text(CP_MJD_arr_cal[idxst][j], ymax * 1.05, \
                                                    target_names_CP_cal[idxst][j].replace('_', ' ')+ \
                                                    " ("+np.str(np.round(target_fluxes_CP_cal[idxst][j],1))+"Jy)", rotation=90,
                                                    va='bottom',fontsize=7)  
                                text_tag_flag = 0
                                prev_text_MJD = CP_MJD_arr_cal[idxst][j]
                                prev_target_name = target_names_CP_cal[idxst][j]

            if len(CP_sta_index) > 0:
                idxst = np.all(CP_sta_index == CP_sta_indices[i], axis=1)
                if len(CP_arr[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(CP_BCD_arr == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(CP_arr[cidxst]) > 0:
                            if plot_errorbars == True:
                                axs[i + 0].errorbar(CP_MJD_arr[cidxst], CP_arr[cidxst], yerr=CPerr_arr[cidxst],
                                                fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                label='CP sci ' + BCD_labels[j])
                            else:
                                axs[i + 0].errorbar(CP_MJD_arr[cidxst], CP_arr[cidxst],
                                                    fmt=BCD_markers[j], color=V2_colors[j], elinewidth=1.5,
                                                    label='CP sci ' + BCD_labels[j])
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if CP_MJD_arr[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_CP[idxst][j]):
                                ymin, ymax = axs[i + 0].get_ylim()
                                axs[i + 0].text(CP_MJD_arr[idxst][j], ymax * 1.05,
                                                target_names_CP[idxst][j].replace('_', ' '),
                                                rotation=90,
                                                va='bottom')
                                text_tag_flag = 0
                                prev_text_MJD = CP_MJD_arr[idxst][j]
                                prev_target_name = target_names_CP[idxst][j]
            axlabel = sta_names[CP_sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                    sta_names[CP_sta_indices[i, 1] == dic['STA_INDEX']][0] + ' - ' + \
                    sta_names[CP_sta_indices[i, 2] == dic['STA_INDEX']][0]
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            axs[i + 0].text(0.05, 0.95, axlabel, horizontalalignment='left',
                            verticalalignment='top',
                            transform=axs[i + 0].transAxes, bbox=props)
            if i == 0:
                leg = axs[i + 0].legend(loc='upper right')
                leg.get_frame().set_alpha(0.5)
            axs[i + 0].set_ylabel(r'$CP\,\left(^\circ\\right)$')
            axs[i + 0].set_xlabel(r'$\mathrm{MJD}$')
        plt.suptitle(r'$CP\mathrm{\ vs.\ time}$')

        fig2.subplots_adjust(hspace=0, wspace=0)
        for i in range(2):
            plt.setp(axs[i + 0].get_xticklabels(), visible=False)
            x_axis = axs[i + 0].axes.get_xaxis()
            x_axis.get_label().set_visible(False)
        for i in range(1, 4, 2):
            plt.setp(axs[i + 0].get_yticklabels(), visible=False)
            y_axis = axs[i + 0].axes.get_yaxis()
            y_axis.get_label().set_visible(False)

        if saveplots:
            fig2.savefig(output_path / f'_CP.png', dpi=150)
            fig2.savefig(output_path / f'_CP.eps', format='eps', dpi=300)
            plt.close(fig2)
        else:
            plt.show()
        print("Plots READY.")


def show_vis_tf_vs_time(
        list_of_dicts: List[Dict],
        wlenRange: List[float],
        saveplots: Optional[bool] = False,
        output_path: Optional[Path] = None,
        plot_errorbars: Optional[bool] = True, **kwargs) -> None:
    """Plots the TF vs. time for the oifits-files.

    Parameters
    ----------
    list_of_dicts : list of dict
        List of oifits-files.
    wlenRange : list
        Range of wavelength to plot.
    saveplots : bool, optional
        If True, the plots are saved in the output_path. The default is False.
    output_path : pathlib.Path, optional
        Path to save the plots. The default is None.
    plot_errorbars : bool, optional
        If True, the errorbars are plotted. The default is True.
    """
    print('Go !')

    if list_of_dicts:
        # colors: BCD: out-out, in-in, in-out, out-in
        BCD_configs = np.array([[0, 0], [1, 1], [1, 0], [0, 1]])
        BCD_labels = np.array(['OUT-OUT', 'IN-IN', 'IN-OUT', 'OUT-IN'])
        BCD_markers = np.array(['o', 's', 'd', 'p'])
        V_cal_colors = np.array(['darkseagreen', 'yellowgreen', 'olivedrab', 'darkkhaki'])
        V_colors = np.array(['red', 'orange', 'salmon', 'pink'])
        TF_colors = np.array(['blue', 'dodgerblue', 'royalblue', 'indigo'])

        target_names_cal = []
        target_fluxes_cal = []
        V_BCD_arr_cal = []
        V_MJD_arr_cal = []
        V_arr_cal = []
        Verr_arr_cal = []
        V_sta_index_cal = []

        target_names_TF = []
        target_fluxes_TF = []
        TF_BCD_arr = []
        TF_MJD_arr = []
        TF_arr = []
        TFerr_arr = []
        TF_sta_index = []

        target_names = []
        target_fluxes = []
        V_BCD_arr = []
        V_MJD_arr = []
        V_arr = []
        Verr_arr = []
        V_sta_index = []

        for dic in list_of_dicts:
            print(dic['TARGET'])
            wl = np.array(dic['WLEN'])
            wlenRange_idx = np.logical_and(wl > wlenRange[0] / 1.0e6, wl < wlenRange[1] / 1.0e6)
            if sum(wlenRange_idx) > 0:
                category = dic['CATEGORY'].lower()
                if 'cal' in category:
                    try:
                        datay = np.sqrt(3.4)*np.array(dic['VIS']['VIS'])
                        datayerr = np.sqrt(3.4)*np.array(dic['VIS']['VISERR'])
                        datax = np.array(dic['VIS']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            V_arr_cal.append(robust.mean(datay[i, wlenRange_idx]))
                            Verr_arr_cal.append(robust.mean(datayerr[i, wlenRange_idx]))
                            V_BCD_arr_cal.append([BCD1, BCD2])
                            V_MJD_arr_cal.append(datax[i])
                            target_names_cal.append(dic['TARGET'])
                            target_fluxes_cal.append(dic['TARGET_FLUX'])
                            sta_index = np.sort(dic['VIS2']['STA_INDEX'][i])
                            V_sta_index_cal.append(sta_index)
                    except Exception:
                        print(dic['TARGET'], dic['DATEOBS'], "No CAL VIS data found.")
                    try:
                        datay = np.sqrt(3.4)*np.array(dic['TF2']['TF'])
                        datayerr = np.sqrt(3.4)*np.array(dic['TF2']['TFERR'])
                        datax = np.array(dic['TF2']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            TF_arr.append(robust.mean(datay[i, wlenRange_idx]))
                            TFerr_arr.append(robust.mean(datayerr[i, wlenRange_idx]))
                            TF_BCD_arr.append([BCD1, BCD2])
                            TF_MJD_arr.append(datax[i])
                            target_names_TF.append(dic['TARGET'])
                            target_fluxes_TF.append(dic['TARGET_FLUX'])
                            sta_index = np.sort(dic['TF2']['STA_INDEX'][i])
                            TF_sta_index.append(sta_index)
                    except Exception:
                        print (dic['TARGET'], dic['DATEOBS'], "No CAL TF data found.")
                if 'sci' in category:
                    try:
                        datay = np.sqrt(3.4)*np.array(dic['VIS']['VIS'])
                        datayerr = np.sqrt(3.4)*np.array(dic['VIS']['VISERR'])
                        datax = np.array(dic['VIS']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            V_arr.append(robust.mean(datay[i, wlenRange_idx]))
                            Verr_arr.append(robust.mean(datayerr[i, wlenRange_idx]))
                            V_BCD_arr.append([BCD1, BCD2])
                            V_MJD_arr.append(datax[i])
                            target_names.append(dic['TARGET'])
                            target_fluxes.append(dic['TARGET_FLUX'])
                            sta_index = np.sort(dic['VIS']['STA_INDEX'][i])
                            V_sta_index.append(sta_index)
                    except Exception:
                        print (dic['TARGET'], dic['DATEOBS'], "No SCI VIS data found.")
            else:
                print("Wavelength out of range.")

        sta_names = dic['STA_NAME']
        target_names_cal = np.array(target_names_cal)
        target_fluxes_cal = np.array(target_fluxes_cal)
        V_BCD_arr_cal = np.array(V_BCD_arr_cal)
        V_MJD_arr_cal = np.array(V_MJD_arr_cal)
        V_arr_cal = np.array(V_arr_cal)
        Verr_arr_cal = np.array(Verr_arr_cal)
        V_sta_index_cal = np.array(V_sta_index_cal)

        target_names_TF = np.array(target_names_TF)
        target_fluxes_TF = np.array(target_fluxes_TF)
        TF_BCD_arr = np.array(TF_BCD_arr)
        TF_MJD_arr = np.array(TF_MJD_arr)
        TF_arr = np.array(TF_arr)
        TFerr_arr = np.array(TFerr_arr)
        TF_sta_index = np.array(TF_sta_index)

        target_names = np.array(target_names)
        target_fluxes = np.array(target_fluxes)
        V_BCD_arr = np.array(V_BCD_arr)
        V_MJD_arr = np.array(V_MJD_arr)
        V_arr = np.array(V_arr)
        Verr_arr = np.array(Verr_arr)
        V_sta_index = np.array(V_sta_index)

        if len(V_sta_index_cal) > 0:
            sta_indices = np.unique(V_sta_index_cal, axis=0)
        elif len(V_sta_index) > 0:
            sta_indices = np.unique(V_sta_index, axis=0)
        else:
            print ("Data arrays empty. Quitting.")
            return
        n_max_config = np.nanmax([6, sta_indices.shape[0]])

        if len(V_MJD_arr_cal) > 0 and len(V_MJD_arr) > 0:
            MJD_range = [np.nanmin([np.nanmin(V_MJD_arr_cal), np.nanmin(V_MJD_arr)]),
                         np.nanmax([np.nanmax(V_MJD_arr_cal), np.nanmax(V_MJD_arr)])]
        elif len(V_MJD_arr) > 0:
            MJD_range = [np.nanmin(V_MJD_arr), np.nanmax(V_MJD_arr)]
        elif len(V_MJD_arr_cal) > 0:
            MJD_range = [np.nanmin(V_MJD_arr_cal), np.nanmax(V_MJD_arr_cal)]
        else:
            MJD_range = [0.0, 1.0]
        text_width_MJD = (MJD_range[1] - MJD_range[0]) / 20.0
        fig1, axs1 = plt.subplots(3, 2, figsize=(15, 16), sharex=True, sharey=True)
        axs1 = axs1.ravel()
        text_y = 1.02
        for i in range(n_max_config):
            if len(V_sta_index_cal) > 0:
                idxst = np.all(V_sta_index_cal == sta_indices[i], axis=1)
                if len(V_arr_cal[idxst]) > 0:
                    label = 'V cal '
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(V_BCD_arr_cal == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)

                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if V_MJD_arr_cal[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_cal[idxst][j]):
                                axs1[i].text(V_MJD_arr_cal[idxst][j], text_y, \
                                                 target_names_cal[idxst][j].replace('_', ' ')+ \
                                                 " ("+np.str(np.round(target_fluxes_cal[idxst][j],1))+"Jy)", rotation=90, \
                                                 va='bottom',fontsize=7)
                                text_tag_flag = 0
                                prev_text_MJD = V_MJD_arr_cal[idxst][j]
                                prev_target_name = target_names_cal[idxst][j]

            if len(TF_sta_index) > 0:
                label = 'TF '
                idxst = np.all(TF_sta_index == sta_indices[i], axis=1)
                if len(TF_arr[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(TF_BCD_arr == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(TF_arr[cidxst]) > 0:
                            if plot_errorbars:
                                axs1[i].errorbar(TF_MJD_arr[cidxst],TF_arr[cidxst],
                                                 yerr=TFerr_arr[cidxst],
                                                 fmt=BCD_markers[j], color=TF_colors[j], elinewidth=1.5,
                                                 label=label + BCD_labels[j])
                                z=np.polyfit(TF_MJD_arr[cidxst]-np.min(TF_MJD_arr[cidxst]), TF_arr[cidxst],3)
                                p=np.poly1d(z)
                                x=(np.max(TF_MJD_arr[cidxst])-np.min(TF_MJD_arr[cidxst]))*np.arange(100)/100.
                                axs1[i].plot(np.min(TF_MJD_arr[cidxst])+x,p(x),color=TF_colors[j])
                                val=np.sqrt(np.abs(TF_arr[cidxst]))-p(TF_MJD_arr[cidxst]-np.min(TF_MJD_arr[cidxst]))
                                valmed=np.median(val)
                                print('MAD={0},{1},{2}'.format(np.median(np.abs(val-valmed)),BCD_labels[j],sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]))
                                print('RMS={0},{1},{2}'.format(np.std(val),BCD_labels[j],sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]))
                                print('Direct RMS={0},{1},{2}'.format(np.std(TF_arr[cidxst]),BCD_labels[j],sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]))
                            else:
                                axs1[i].errorbar(TF_MJD_arr[cidxst], TF_arr[cidxst],
                                                 fmt=BCD_markers[j], color=TF_colors[j], elinewidth=1.5,
                                                 label=label + BCD_labels[j])
            if len(V_sta_index) > 0:
                label = 'V sci '
                idxst = np.all(V_sta_index == sta_indices[i], axis=1)
                if len(V_arr[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(V_BCD_arr == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(V_arr[cidxst]) > 0:
                            if plot_errorbars:
                                axs1[i].errorbar(V_MJD_arr[cidxst], V_arr[cidxst], yerr=Verr_arr[cidxst],
                                                 fmt=BCD_markers[j], color=V_colors[j], elinewidth=1.5,
                                                 label=label + BCD_labels[j])
                            else:
                                axs1[i].errorbar(V_MJD_arr[cidxst], V_arr[cidxst],
                                                 fmt=BCD_markers[j], color=V_colors[j], elinewidth=1.5,
                                                 label=label + BCD_labels[j])
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if V_MJD_arr[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names[idxst][j]):
                                axs1[i].text(V_MJD_arr[idxst][j], text_y, target_names[idxst][j].replace('_', ' '),
                                             rotation=90, va='bottom', color='darkred')
                                text_tag_flag = 0
                                prev_text_MJD = V_MJD_arr[idxst][j]
                                prev_target_name = target_names[idxst][j]

            axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                    sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            axs1[i].text(0.05, 0.95, axlabel, horizontalalignment='left', verticalalignment='top',
                         transform=axs1[i].transAxes, bbox=props)
            if i == 0:
                leg = axs1[i].legend(loc='upper right')
                leg.get_frame().set_alpha(0.5)
            axs1[i].set_ylim([0., 1.0])
            ylabel = '$V$'
            axs1[i].set_ylabel(ylabel)
            axs1[i].set_xlabel(r'$\mathrm{MJD}$')

        plt.suptitle(r'$V\mathrm{\ vs.\ time}$')
        fig1.subplots_adjust(hspace=0, wspace=0)
        for i in range(4):
            plt.setp(axs1[i].get_xticklabels(), visible=False)
            x_axis = axs1[i].axes.get_xaxis()
            x_axis.get_label().set_visible(False)
        for i in range(1, 6, 2):
            plt.setp(axs1[i].get_yticklabels(), visible=False)
            y_axis = axs1[i].axes.get_yaxis()
            y_axis.get_label().set_visible(False)

        if saveplots:
            label = '_VIS_TF'
            fig1.savefig(output_path / f'{label}.png', dpi=150)
            fig1.savefig(output_path / f'{label}.eps', format='eps', dpi=300)
            plt.close(fig1)
        else:
            plt.show()
        print ("Plots READY.")
        

def show_cf2_vs_time(
        list_of_dicts: List[Dict],
        wlenRange: List[float],
        showvis: Optional[bool] = False,
        saveplots: Optional[bool] = False,
        output_path: Optional[Path] = None, 
        plot_errorbars: Optional[bool] = True):
    """Plots the cf2 vs time.

    Parameters
    ----------
    list_of_dicts : list of dict
        List of oifits-files.
    wlenRange : list
        Wavelength range for the plot.
    showvis : bool, optional
        If True, plot visibilities (V) instead of V^2: V is calculated from V^2 (not from the VISAMP table). The default is False.
    saveplots : bool, optional
        If True, the plots are saved in the output_path. The default is False.
    output_path : pathlib.Path, optional
        Path to save the plots. The default is None.
    plot_errorbars : bool, optional
        If True, plot errorbars. The default is True.
    """
    #starFlux={'Teta Pyx':27, 'psi Vir':44,'tet Cen':48, \
    #      'del Vir':160,'eps Sco':52,'del Oph':110, \
    #      'FZ Lib':17,'del Sco':15,'HD171094':34, \
    #      'e Aql':27,'29 Cap':31,'alf Ant':19, \
    #      'HD99333':20,'C01 Cen':22,'alpha Arae':10, \
    #      'HD138505':18,'HD142527':12,'HD161560':0.4, \
    #      'HD171094':34,'AV Mic':18,'HD181041':6,  \
    #      'HD138492':8,'HD147929':8,'HD148255':3, \
    #      'G Sco':26,'HD156936':3,'HD161849':8, \
    #      'HD165413':5,'HD186765':9,'Mu Hya':30, \
    #      'CH Vir':20,'HD 126111':6,'LY TrA':7, \
    #      'ET Vir':29,'H Sco':30,'HD177507':17, \
    #      'V345 Tel':9,'RW Lup':18,'HD 150798':144, \
    #      'BM Sco':85,'RX Tel':82,'HD328913':100, \
    #      'Eps Oph':16,'HD182669':2.5,'Nu Hya':30}    # check if list is not empty:
    # dicoFluxL={'lamPyx': 28.105587, 'HD56618': 196.27484, 'omeLup': 115.3, 'Khya': 60.099998, 'AlfCMa': 995.79999, 'HD92397': 129.69469, 'epsCru': 217.05141, 'HD117597': 8.8213739, 'HD151843': 2.9975743, 'etaCen': 20.529898, 'delCen': 30.857565, 'HD100546': 20.5, 'HD107625': 0.49486896, 'BCen': 72.300003, 'HD162415': 12.079382, 'HD114873': 20.785559, 'HD103515': 24.602552, 'HD104329': 2.9723437, 'HD107869': 29.714039, 'muCen': 9.1099892, 'cVel': 108.7, 'tetCMa': 164.47496, 'HD109960': 17.799999, 'HD148984': 14.406607, 'HD150071': 0.31112576, 'HD105381': 12.466991, 'delSco': 32.341614, 'V1068Sco': 317.78168, 'HSco': 191.98087, 'HD140143': 9.6999998, 'HD97705': 0.043862239, 'HD148441': 0.99144912, '27Sco': 97.5, 'TetCen': 377.60001, 'HD102766': 145.40553, 'HD103975': 2.9723437, 'HD108110': 13.865284, 'HD97765': 0.12568563, 'HD113601': 5.0366859, 'HD148730': 0.20337902, '3C273': 0.059742939,'cVel': 108.7, 'omeLup': 115.3, 'Khya': 60.099998, 'AlfCMa': 995.79999, 'HD92397': 129.69469, 'epsCru': 217.05141, 'delCen': 30.857565, 'HSco': 191.98087, 'HD94495': 30.6, 'HD150071': 0.31112576, 'HD162189': 75.758179, 'BCen': 72.300003, 'TetCen': 377.60001, 'HD100546': 20.5, 'HD103515': 24.602552, 'HD97765': 0.12568563}
    # dicoFluxN= {'lamPyx': 4.1238999, 'HD56618': 30.814779, 'omeLup': 17.6, 'Khya': 8.5500002, 'AlfCMa': 118.0, 'HD92397': 20.248549, 'epsCru': 32.011219, 'HD117597': 1.3875213, 'HD151843': 0.51499999, 'etaCen': 5.8867059, 'delCen': 14.54, 'HD100546': 45.694008, 'HD107625': 0.054777335, 'BCen': 10.616354, 'HD162415': 1.9299999, 'HD114873': 3.0, 'HD103515': 5.0599999, 'HD104329': 0.42240164, 'HD107869': 5.4688506, 'muCen': 2.6417792, 'cVel': 15.4, 'tetCMa': 30.499096, 'HD109960': 2.9200001, 'HD148984': 2.2, 'HD150071': 0.032316711, 'HD105381': 1.6773663, 'delSco': 12.018785, 'V1068Sco': 50.355038, 'HSco': 30.69574, 'HD140143': 1.71, 'HD97705': 0.0048251506, 'HD148441': 0.10740463, '27Sco': 10.9, 'TetCen': 48.400002, 'HD102766': 37.23, 'HD103975': 0.42240164, 'HD108110': 2.3499999, 'HD97765': 0.013920446, 'HD113601': 0.76300001, 'HD148730': 0.022324009, '3C273': 0.0067190719,'cVel': 15.4, 'omeLup': 17.6, 'Khya': 8.5500002, 'AlfCMa': 118.0, 'HD92397': 20.248549, 'epsCru': 32.011219, 'delCen': 14.54, 'HSco': 30.69574, 'HD94495': 5.02, 'HD150071': 0.032316711, 'HD162189': 13.1, 'BCen': 10.616354, 'TetCen': 48.400002, 'HD100546': 45.694008, 'HD103515': 5.0599999, 'HD97765': 0.013920446}
    starFlux={'lamPyx': 16.899242, 'HD56618': 123.59395, 'omeLup': 64.288948, 'Khya': 34.413467, 'AlfCMa': 535.90002, 'HD92397': 79.279999, 'epsCru': 143.89999, 'HD117597': 5.3375411, 'HD151843': 1.8043216, 'etaCen': 18.671946, 'delCen': 27.200001, 'HD100546': 13.9, 'HD107625': 0.2566978, 'BCen': 40.465576, 'HD162415': 10.16, 'HD114873': 13.8, 'HD103515': 23.01, 'HD104329': 1.7560067, 'HD107869': 17.4, 'muCen': 6.5473495, 'cVel': 67.809998, 'tetCMa': 92.122192, 'HD109960': 9.8630896, 'HD148984': 5.6999998, 'HD150071': 0.16719383, 'HD105381': 5.1822686, 'delSco': 32.172329, 'V1068Sco': 198.0, 'HSco': 101.3, 'HD140143': 3.5664213, 'HD97705': 0.024131674, 'HD148441': 0.55915225, '27Sco': 48.699944, 'TetCen': 184.3, 'HD102766': 52.769211, 'HD103975': 1.7560067, 'HD108110': 8.4795961, 'HD97765': 0.068838015, 'HD113601': 4.2104497, 'HD148730': 0.11200868, '3C273': 0.03288243,'cVel': 67.809998, 'omeLup': 64.288948, 'Khya': 34.413467, 'AlfCMa': 535.90002, 'HD92397': 79.279999, 'epsCru': 143.89999, 'delCen': 27.200001, 'HSco': 101.3, 'HD94495': 12.719484, 'HD150071': 0.16719383, 'HD162189': 37.186783, 'BCen': 40.465576, 'TetCen': 184.3, 'HD100546': 13.9, 'HD103515': 23.01, 'HD97765': 0.068838015}

    # ('L', {'omeLup': 115.3, 'HD56618': 196.27484, 'HD148730': 0.20337902, 'Khya': 60.099998, 'AlfCMa': 995.79999, 'HD148984': 14.406607, 'HSco': 191.98087, 'epsCru': 217.05141, 'HD148441': 0.99144912, 'delCen': 30.857565, 'HD107869': 29.714039, 'HD104329': 2.9723437, 'HD105381': 12.466991, 'HD107625': 0.49486896, 'BCen': 72.300003, 'TetCen': 377.60001, 'HD103515': 24.602552, 'HD108110': 13.865284, 'HD140143': 9.6999998})
    #('N', {'omeLup': 17.6, 'HD56618': 30.814779, 'HD148730': 0.022324009, 'Khya': 8.5500002, 'AlfCMa': 118.0, 'HD148984': 2.2, 'HSco': 30.69574, 'epsCru': 32.011219, 'HD148441': 0.10740463, 'delCen': 14.54, 'HD107869': 5.4688506, 'HD104329': 0.42240164, 'HD105381': 1.6773663, 'HD107625': 0.054777335, 'BCen': 10.616354, 'TetCen': 48.400002, 'HD103515': 5.0599999, 'HD108110': 2.3499999, 'HD140143': 1.71})
    #('M', {'omeLup': 64.288948, 'HD56618': 123.59395, 'HD148730': 0.11200868, 'Khya': 34.413467, 'AlfCMa': 535.90002, 'HD148984': 5.6999998, 'HSco': 101.3, 'epsCru': 143.89999, 'HD148441': 0.55915225, 'delCen': 27.200001, 'HD107869': 17.4, 'HD104329': 1.7560067, 'HD105381': 5.1822686, 'HD107625': 0.2566978, 'BCen': 40.465576, 'TetCen': 184.3, 'HD103515': 23.01, 'HD108110': 8.4795961, 'HD140143': 3.5664213})

    if list_of_dicts:
        # colors: BCD: out-out, in-in, in-out, out-in
        BCD_configs = np.array([[0, 0], [1, 1], [1, 0], [0, 1]])
        BCD_labels = np.array(['OUT-OUT', 'IN-IN', 'IN-OUT', 'OUT-IN'])
        BCD_markers = np.array(['o', 's', 'd', 'p'])
        CF2_cal_colors = np.array(['darkseagreen', 'yellowgreen', 'olivedrab', 'darkkhaki'])
        CF2_colors = np.array(['red', 'orange', 'salmon', 'pink'])

        target_names_cal = []
        CF2_BCD_arr_cal = []
        CF2_MJD_arr_cal = []
        CF2_arr_cal = []
        CF2err_arr_cal = []
        CF2_sta_index_cal = []

        target_names_CP_cal = []
        CP_BCD_arr_cal = []
        CP_MJD_arr_cal = []
        CP_arr_cal = []
        CPerr_arr_cal = []
        CP_sta_index_cal = []

        target_names = []
        CF2_BCD_arr = []
        CF2_MJD_arr = []
        CF2_arr = []
        CF2err_arr = []
        CF2_sta_index = []

        target_names_CP = []
        CP_BCD_arr = []
        CP_MJD_arr = []
        CP_arr = []
        CPerr_arr = []
        CP_sta_index = []

        for dic in list_of_dicts:
            print(dic['TARGET'])
            wl = np.array(dic['WLEN'])
            wlenRange_idx = np.logical_and(wl > wlenRange[0] / 1.0e6, wl < wlenRange[1] / 1.0e6)
            if sum(wlenRange_idx) > 0:
                category = dic['CATEGORY'].lower()
                if 'cal' in category:
                    try:
                        datay = np.array(dic['CF2']['CF2_RATIO'])
                        datayerr = np.array(dic['CF2']['CF2_RATIOERR'])
                        datax = np.array(dic['CF2']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            CF2_arr_cal.append(robust.mean(datay[i, wlenRange_idx]))
                            CF2err_arr_cal.append(robust.mean(datayerr[i, wlenRange_idx]))
                            CF2_BCD_arr_cal.append([BCD1, BCD2])
                            CF2_MJD_arr_cal.append(datax[i])
                            target_names_cal.append(dic['TARGET'])
                            sta_index = np.sort(dic['CF2']['STA_INDEX'][i])
                            CF2_sta_index_cal.append(sta_index)
                    except Exception:
                        print (dic['TARGET'], dic['DATEOBS'], "No CAL CF2 ratio data found.")
                    try:
                        datay = np.array(dic['T3']['CLOS'])
                        datayerr = np.array(dic['T3']['CLOSERR'])
                        datax = np.array(dic['T3']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            CP_arr_cal.append(robust.mean(datay[i, wlenRange_idx]))
                            CPerr_arr_cal.append(robust.mean(datayerr[i, wlenRange_idx]))
                            CP_BCD_arr_cal.append([BCD1, BCD2])
                            CP_MJD_arr_cal.append(datax[i])
                            target_names_CP_cal.append(dic['TARGET'])
                            sta_index = np.sort(dic['T3']['STA_INDEX'][i])
                            CP_sta_index_cal.append(sta_index)
                    except Exception:
                        print (dic['TARGET'], dic['DATEOBS'], "No CAL CP data found.")
                if 'sci' in category:
                    try:
                        datay = np.array(dic['CF2']['CF2_RATIO'])
                        datayerr = np.array(dic['CF2']['CF2_RATIOERR'])
                        datax = np.array(dic['CF2']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            CF2_arr.append(robust.mean(datay[i, wlenRange_idx]))
                            CF2err_arr.append(robust.mean(datayerr[i, wlenRange_idx]))
                            CF2_BCD_arr.append([BCD1, BCD2])
                            CF2_MJD_arr.append(datax[i])
                            target_names.append(dic['TARGET'])
                            sta_index = np.sort(dic['CF2']['STA_INDEX'][i])
                            CF2_sta_index.append(sta_index)
                    except:
                        print (dic['TARGET'], dic['DATEOBS'], "No SCI CF2 ratio data found.")
                    try:
                        datay = np.array(dic['T3']['CLOS'])
                        datayerr = np.array(dic['T3']['CLOSERR'])
                        datax = np.array(dic['T3']["TIME"])
                        n_rows = datay.shape[0]
                        for i in range(n_rows):
                            if dic['BCD1NAME'] == 'IN':
                                BCD1 = 1
                            elif dic['BCD1NAME'] == 'OUT':
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic['BCD2NAME'] == 'IN':
                                BCD2 = 1
                            elif dic['BCD2NAME'] == 'OUT':
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            CP_arr.append(robust.mean(datay[i, wlenRange_idx]))
                            CPerr_arr.append(robust.mean(datayerr[i, wlenRange_idx]))
                            target_names_CP.append(dic['TARGET'])
                            sta_index = np.sort(dic['T3']['STA_INDEX'][i])
                            CP_sta_index.append(sta_index)
                            CP_BCD_arr.append([BCD1, BCD2])
                            CP_MJD_arr.append(datax[i])
                    except:
                        print (dic['TARGET'], dic['DATEOBS'], "No SCI CP data found.")
            else:
                print("Wavelength out of range.")

        sta_names = dic['STA_NAME']
        target_names_cal = np.array(target_names_cal)
        CF2_BCD_arr_cal = np.array(CF2_BCD_arr_cal)
        CF2_MJD_arr_cal = np.array(CF2_MJD_arr_cal)
        CF2_arr_cal = np.array(CF2_arr_cal)
        CF2err_arr_cal = np.array(CF2err_arr_cal)
        CF2_sta_index_cal = np.array(CF2_sta_index_cal)

        target_names_CP_cal = np.array(target_names_CP_cal)
        CP_BCD_arr_cal = np.array(CP_BCD_arr_cal)
        CP_MJD_arr_cal = np.array(CP_MJD_arr_cal)
        CP_arr_cal = np.array(CP_arr_cal)
        CPerr_arr_cal = np.array(CPerr_arr_cal)
        CP_sta_index_cal = np.array(CP_sta_index_cal)

        target_names = np.array(target_names)
        CF2_BCD_arr = np.array(CF2_BCD_arr)
        CF2_MJD_arr = np.array(CF2_MJD_arr)
        CF2_arr = np.array(CF2_arr)
        CF2err_arr = np.array(CF2err_arr)
        CF2_sta_index = np.array(CF2_sta_index)

        target_names_CP = np.array(target_names_CP)
        CP_BCD_arr = np.array(CP_BCD_arr)
        CP_MJD_arr = np.array(CP_MJD_arr)
        CP_arr = np.array(CP_arr)
        CPerr_arr = np.array(CPerr_arr)
        CP_sta_index = np.array(CP_sta_index)

        if len(CF2_sta_index_cal) > 0:
            sta_indices = np.unique(CF2_sta_index_cal, axis=0)
        elif len(CF2_sta_index) > 0:
            sta_indices = np.unique(CF2_sta_index, axis=0)
        else:
            print ("Data arrays empty. Quitting.")
            return
        n_max_config = np.nanmax([6, sta_indices.shape[0]])
        print('sta_indices[:,0] = {0}'.format(sta_indices[:,0]))
        print('sta_indices[:,1] = {0}'.format(sta_indices[:,1]))
        print(dic['STA_INDEX'])

        if len(CF2_MJD_arr_cal) > 0 and len(CF2_MJD_arr) > 0:
            MJD_range = [np.nanmin([np.nanmin(CF2_MJD_arr_cal), np.nanmin(CF2_MJD_arr)]),
                         np.nanmax([np.nanmax(CF2_MJD_arr_cal), np.nanmax(CF2_MJD_arr)])]
        elif len(CF2_MJD_arr) > 0:
            MJD_range = [np.nanmin(CF2_MJD_arr), np.nanmax(CF2_MJD_arr)]
        elif len(CF2_MJD_arr_cal) > 0:
            MJD_range = [np.nanmin(CF2_MJD_arr_cal), np.nanmax(CF2_MJD_arr_cal)]
        else:
            MJD_range = [0.0, 1.0]
        text_width_MJD = (MJD_range[1] - MJD_range[0]) / 20.0
        fig1, axs1 = plt.subplots(3, 2, figsize=(15, 16), sharex=True, sharey=True)
        axs1 = axs1.ravel()
        text_y = 2.05
        for i in range(n_max_config):
            if len(CF2_sta_index_cal) > 0:
                idxst = np.all(CF2_sta_index_cal == sta_indices[i], axis=1)
                if len(CF2_arr_cal[idxst]) > 0:
                    if showvis:
                        label = 'CF_ratio cal '
                    else:
                        label = 'CF2_ratio cal '
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(CF2_BCD_arr_cal == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        print('i = {0}'.format(i))
                        if len(CF2_arr_cal[cidxst]) > 0:
                            if plot_errorbars:
                                if showvis:
                                    axs1[i].errorbar(CF2_MJD_arr_cal[cidxst], np.sqrt(np.abs(CF2_arr_cal[cidxst])),
                                                     yerr=0.5 * CF2err_arr_cal[cidxst] / np.sqrt(np.abs(CF2_arr_cal[cidxst])),
                                                     fmt=BCD_markers[j], color=CF2_cal_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                    z=np.polyfit(CF2_MJD_arr_cal[cidxst]-np.min(CF2_MJD_arr_cal[cidxst]), np.sqrt(np.abs(CF2_arr_cal[cidxst])),3)
                                    p=np.poly1d(z)
                                    x=(np.max(CF2_MJD_arr_cal[cidxst])-np.min(CF2_MJD_arr_cal[cidxst]))*np.arange(100)/100.
                                    axs1[i].plot(np.min(CF2_MJD_arr_cal[cidxst])+x,p(x),color=CF2_cal_colors[j])
                                    val=np.sqrt(np.abs(CF2_arr_cal[cidxst]))-p(CF2_MJD_arr_cal[cidxst]-np.min(CF2_MJD_arr_cal[cidxst]))
                                    valmed=np.median(val)

                                    print('MAD={0},{1},{2}'.format(np.median(np.abs(val-valmed)),BCD_labels[j],sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]))
                                    print('RMS={0},{1},{2}'.format(np.std(val),BCD_labels[j],sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]))
                                    
                                else:
                                    axs1[i].errorbar(CF2_MJD_arr_cal[cidxst], CF2_arr_cal[cidxst],
                                                     yerr=CF2err_arr_cal[cidxst],
                                                     fmt=BCD_markers[j], color=CF2_cal_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                            else:
                                if showvis:
                                    axs1[i].errorbar(CF2_MJD_arr_cal[cidxst], np.sqrt(np.abs(CF2_arr_cal[cidxst])),
                                                     fmt=BCD_markers[j], color=CF2_cal_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                else:
                                    axs1[i].errorbar(CF2_MJD_arr_cal[cidxst], CF2_arr_cal[cidxst],
                                                     fmt=BCD_markers[j], color=CF2_cal_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if CF2_MJD_arr_cal[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_cal[idxst][j]):
                                tar_name=target_names_cal[idxst][j].replace('_','').replace(' ','')
                                if (tar_name in starFlux) :
                                    axs1[i].text(CF2_MJD_arr_cal[idxst][j], text_y, \
                                                 target_names_cal[idxst][j].replace('_', ' ')+ \
                                                 " ("+np.str(np.round(starFlux[tar_name],1))+"Jy)", rotation=90, \
                                                 va='bottom',fontsize=7)
                                else :
                                    axs1[i].text(CF2_MJD_arr_cal[idxst][j], text_y, \
                                                 target_names_cal[idxst][j].replace('_', ' '), rotation=90, \
                                                 va='bottom',fontsize=7)
                                    
 
                                text_tag_flag = 0
                                prev_text_MJD = CF2_MJD_arr_cal[idxst][j]
                                prev_target_name = target_names_cal[idxst][j]

            if len(CF2_sta_index) > 0:
                if showvis:
                    label = 'CF ratio sci '
                else:
                    label = 'CF2 ratio sci '
                idxst = np.all(CF2_sta_index == sta_indices[i], axis=1)
                if len(CF2_arr[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(CF2_BCD_arr == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(V2_arr[cidxst]) > 0:
                            if plot_errorbars:
                                if showvis:
                                    axs1[i].errorbar(CF2_MJD_arr[cidxst], np.sqrt(CF2_arr[cidxst]),
                                                     yerr=0.5 * CF2err_arr[cidxst] / np.sqrt(CF2_arr[cidxst]),
                                                     fmt=BCD_markers[j], color=CF2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                else:
                                    axs1[i].errorbar(CF2_MJD_arr[cidxst], CF2_arr[cidxst], yerr=CF2err_arr[cidxst],
                                                     fmt=BCD_markers[j], color=CF2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                            else:
                                if showvis:
                                    axs1[i].errorbar(V2_MJD_arr[cidxst], np.sqrt(CF2_arr[cidxst]),
                                                     fmt=BCD_markers[j], color=CF2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                                else:
                                    axs1[i].errorbar(CF2_MJD_arr[cidxst], CF2_arr[cidxst],
                                                     fmt=BCD_markers[j], color=CF2_colors[j], elinewidth=1.5,
                                                     label=label + BCD_labels[j])
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if CF2_MJD_arr[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names[idxst][j]):
                                axs1[i].text(CF2_MJD_arr[idxst][j], text_y, target_names[idxst][j].replace('_', ' '),
                                             rotation=90, va='bottom', color='darkred')
                                text_tag_flag = 0
                                prev_text_MJD = CF2_MJD_arr[idxst][j]
                                prev_target_name = target_names[idxst][j]

            axlabel = sta_names[sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                      sta_names[sta_indices[i, 1] == dic['STA_INDEX']][0]
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            axs1[i].text(0.05, 0.95, axlabel, horizontalalignment='left', verticalalignment='top',
                         transform=axs1[i].transAxes, bbox=props)
            if i == 0:
                leg = axs1[i].legend(loc='upper right')
                leg.get_frame().set_alpha(0.5)
            axs1[i].set_ylim([-0.1, 2.0])
            if showvis:
                ylabel = '$CF ratio$'
            else:
                ylabel = '$CF^2 ratio$'
            axs1[i].set_ylabel(ylabel)
            axs1[i].set_xlabel(r'$\mathrm{MJD}$')
        if showvis:
            plt.suptitle(r'$CF ratio\mathrm{\ vs.\ time}$')
        else:
            plt.suptitle(r'$CF^2 ratio\mathrm{\ vs.\ time}$')
        fig1.subplots_adjust(hspace=0, wspace=0)

        for i in range(4):
            plt.setp(axs1[i].get_xticklabels(), visible=False)
            x_axis = axs1[i].axes.get_xaxis()
            x_axis.get_label().set_visible(False)

        for i in range(1, 6, 2):
            plt.setp(axs1[i].get_yticklabels(), visible=False)
            y_axis = axs1[i].axes.get_yaxis()
            y_axis.get_label().set_visible(False)
        if saveplots == True:
            if showvis == True:
                label = '_CF_RATIO'
            else:
                label = '_CF2_RATIO'
            fig1.savefig(output_path + label + '.png', dpi=150)
            fig1.savefig(output_path + label + '.eps', format='eps', dpi=300)
            plt.close(fig1)
            
        if len(CP_sta_index_cal) > 0:
            CP_sta_indices = np.unique(CP_sta_index_cal, axis=0)
        else:
            CP_sta_indices = np.unique(CP_sta_index, axis=0)

        # print CP_sta_indices
        n_max_config = np.nanmax([4, CP_sta_indices.shape[0]])

        fig2, axs = plt.subplots(2, 2, figsize=(14, 12), sharex=True, sharey=True)
        axs = axs.ravel()
        text_y = 60

        for i in range(n_max_config):
            axs[i + 0].plot(MJD_range, [0.0, 0.0], '-', color='gray', lw=1.5)
            if len(CP_sta_index_cal) > 0:
                idxst = np.all(CP_sta_index_cal == CP_sta_indices[i], axis=1)
                print(np.shape(CP_MJD_arr_cal))
                print(np.shape(CP_arr_cal))
                if len(CP_arr_cal[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(CP_BCD_arr_cal == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(CP_arr_cal[cidxst]) > 0:
                            if plot_errorbars == True:
                                axs[i + 0].errorbar(CP_MJD_arr_cal[cidxst], CP_arr_cal[cidxst], yerr=CPerr_arr_cal[cidxst],
                                                fmt=BCD_markers[j], color=CF2_cal_colors[j], elinewidth=1.5,
                                                label='CP cal ' + BCD_labels[j])
                                axs[i+0].set_ylim(-10,10)
                                z=np.polyfit(CP_MJD_arr_cal[cidxst]-np.min(CP_MJD_arr_cal[cidxst]), CP_arr_cal[cidxst],3)
                                p=np.poly1d(z)
                                #print('p = {0}').format(p)
                                x=(np.max(CP_MJD_arr_cal[cidxst])-np.min(CP_MJD_arr_cal[cidxst]))*np.arange(100)/100.
                                axs[i].plot(np.min(CP_MJD_arr_cal[cidxst])+x,p(x),color=CF2_cal_colors[j])
                                val=CP_arr_cal[cidxst]-p(CP_MJD_arr_cal[cidxst]-np.min(CP_MJD_arr_cal[cidxst]))
                                valmed=np.median(val)
                                print('MAD={0},{1},{2}'.format(np.median(np.abs(val-valmed)),BCD_labels[j],sta_names[CP_sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[CP_sta_indices[i, 1] == dic['STA_INDEX']][0]+ ' - ' + sta_names[CP_sta_indices[i, 2] == dic['STA_INDEX']][0]))
                                print('RMS={0},{1},{2}'.format(np.std(val),BCD_labels[j],sta_names[CP_sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + sta_names[CP_sta_indices[i, 1] == dic['STA_INDEX']][0]+ ' - ' + sta_names[CP_sta_indices[i, 2] == dic['STA_INDEX']][0]))
                            else:
                                axs[i + 0].errorbar(CP_MJD_arr_cal[cidxst], CP_arr_cal[cidxst],
                                                    fmt=BCD_markers[j], color=CF2_cal_colors[j], elinewidth=1.5,
                                                    label='CP cal ' + BCD_labels[j])
                                axs[i+0].set_ylim(-20,20)
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if CP_MJD_arr_cal[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_CP_cal[idxst][j]):
                                ymin, ymax = axs[i + 0].get_ylim()
                                tar_name=target_names_CP_cal[idxst][j].replace('_','').replace(' ','')
                                if (tar_name in starFlux):
                                    axs[i + 0].text(CP_MJD_arr_cal[idxst][j], ymax * 1.05, \
                                                    target_names_CP_cal[idxst][j].replace('_', ' ')+ \
                                                    " ("+np.str(np.round(starFlux[tar_name],1))+"Jy)", rotation=90,
                                                    va='bottom',fontsize=7)
                                else:
                                     axs[i + 0].text(CP_MJD_arr_cal[idxst][j], ymax * 1.05, \
                                                    target_names_CP_cal[idxst][j].replace('_', ' '), rotation=90,
                                                     va='bottom',fontsize=7)
                                   
                                    
                                text_tag_flag = 0
                                prev_text_MJD = CP_MJD_arr_cal[idxst][j]
                                prev_target_name = target_names_CP_cal[idxst][j]
            if len(CP_sta_index) > 0:
                idxst = np.all(CP_sta_index == CP_sta_indices[i], axis=1)
                if len(CP_arr[idxst]) > 0:
                    for j in range(len(BCD_configs)):
                        BCDidx = np.all(CP_BCD_arr == BCD_configs[j], axis=1)
                        cidxst = np.logical_and(idxst, BCDidx)
                        if len(CP_arr[cidxst]) > 0:
                            if plot_errorbars == True:
                                axs[i + 0].errorbar(CP_MJD_arr[cidxst], CP_arr[cidxst], yerr=CPerr_arr[cidxst],
                                                fmt=BCD_markers[j], color=CF2_colors[j], elinewidth=1.5,
                                                label='CP sci ' + BCD_labels[j])
                            else:
                                axs[i + 0].errorbar(CP_MJD_arr[cidxst], CP_arr[cidxst],
                                                    fmt=BCD_markers[j], color=CF2_colors[j], elinewidth=1.5,
                                                    label='CP sci ' + BCD_labels[j])
                    if i in range(2):
                        text_tag_flag = 1
                        prev_text_MJD = 0.0
                        prev_target_name = ""
                        for j in range(np.sum(idxst)):
                            if CP_MJD_arr[idxst][j] > (prev_text_MJD + text_width_MJD):
                                text_tag_flag = 1
                            if text_tag_flag == 1 or (prev_target_name != target_names_CP[idxst][j]):
                                ymin, ymax = axs[i + 0].get_ylim()
                                axs[i + 0].text(CP_MJD_arr[idxst][j], ymax * 1.05,
                                                target_names_CP[idxst][j].replace('_', ' '),
                                                rotation=90,
                                                va='bottom')
                                text_tag_flag = 0
                                prev_text_MJD = CP_MJD_arr[idxst][j]
                                prev_target_name = target_names_CP[idxst][j]
            axlabel = sta_names[CP_sta_indices[i, 0] == dic['STA_INDEX']][0] + ' - ' + \
                      sta_names[CP_sta_indices[i, 1] == dic['STA_INDEX']][0] + ' - ' + \
                      sta_names[CP_sta_indices[i, 2] == dic['STA_INDEX']][0]
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            axs[i + 0].text(0.05, 0.95, axlabel, horizontalalignment='left', verticalalignment='top',
                            transform=axs[i + 0].transAxes, bbox=props)
            if i == 0:
                leg = axs[i + 0].legend(loc='upper right')
                leg.get_frame().set_alpha(0.5)
            axs[i + 0].set_ylabel(r'$CP\,\left(^\circ\\right)$')
            axs[i + 0].set_xlabel(r'$\mathrm{MJD}$')
        plt.suptitle(r'$CP\mathrm{\ vs.\ time}$')

        fig2.subplots_adjust(hspace=0, wspace=0)
        for i in range(2):
            plt.setp(axs[i + 0].get_xticklabels(), visible=False)
            x_axis = axs[i + 0].axes.get_xaxis()
            x_axis.get_label().set_visible(False)
        for i in range(1, 4, 2):
            plt.setp(axs[i + 0].get_yticklabels(), visible=False)
            y_axis = axs[i + 0].axes.get_yaxis()
            y_axis.get_label().set_visible(False)

        if saveplots:
            fig2.savefig(output_path / f'_CP.png', dpi=150)
            fig2.savefig(output_path / f'_CP.eps', format='eps', dpi=300)
            plt.close(fig2)
        else:
            plt.show()
        print("Plots READY.")


def open_oi_dir(input_dir: Path, choice_band_LM: str):
    """Opens all the oifits-files in the directory."""
    oifits_file_list = list(map(str, Path(input_dir).glob('*RAW_INT*fits')))

    list_of_dicts = []
    for file in oifits_file_list:
        if "LAMP" not in file:
            dic = open_oi(file, choice_band_LM)
            if dic:
                print (dic['TARGET'], dic['DATEOBS'], dic['BAND'], dic['DISP'], dic['DIT'], dic['CATEGORY'])
                list_of_dicts.append(dic)
    return list_of_dicts


def filter_oi_list(
        list_of_dicts: List[Dict],
        dates: Optional[List[str]] = [],
        bands: Optional[List[str]] = [],
        spectral_resolutions: Optional[List[str]] = [],
        DIT_range: Optional[List[float]] = [],
        targets: Optional[List[str]] = [],
        bcd_config: Optional[List[str]] = [],
        seeing_range: Optional[List[float]] = [],
        flux_range: Optional[List[float]] = []) -> List[Dict]:
    """Filters the list of oifits-files.

    Parameters
    ----------
    list_of_dicts : list of dict
        List of oifits-files.
    dates : str or list of str
        Date of the observation. The format is ["YYYY-MM-DD"].
        The default is '2000-01-01'.
    bands : list
        List of bands. 'L','M','LM', 'N'. The default is [].
    spectral_resolutions : list of str
        List of spectral resolutions. 'LOW','MED','HIGH'.
        The default is [].
    DIT_range : list of float
        List of DIT-range. The format is [min,max] (s).
        The default is [].
    targets : list of str
        List of targets. The default is [].
    bcd_config : list of str
        List of bcd-config. The default is [].
    seeing_range : list of float
        List of seeing-range. The default is [].
    flux_range : list of float
        List of flux-range. The default is [].

    Returns
    -------
    list of dict
        List of filtered oifits-files.
    """
    filtered_list_of_dicts = []
    if bands:
        bands_new = []
        for i in range(len(bands)):
            if bands[i] == 'M':
                bands_new.append('LM')
            elif bands[i] == 'L':
                bands_new.append('LM')
            else:
                bands_new.append(bands[i])

    seeing_value=[]
    tau0_value=[]
    for dic in list_of_dicts:
        if dic:
            date = dic['DATEOBS'][0:10]
            if dates:
                if date not in dates:
                    continue
            if bands:
                if dic['BAND'] not in bands_new:
                    breakpoint()
                    continue
            if spectral_resolutions:
                if dic['DISP'] not in spectral_resolutions:
                    continue
            if bcd_config:
                if dic['BCDCONFIG'] not in bcd_config:
                    continue
            if DIT_range:
                if not (dic['DIT'] >= DIT_range[0] and dic['DIT'] <= DIT_range[1]):
                    continue
            if seeing_range:
                if not (dic['SEEING'] >= seeing_range[0] and dic['SEEING'] <= seeing_range[1]):
                    continue
            if flux_range:
                if not (dic['TARGET_FLUX'] >= flux_range[0] and dic['TARGET_FLUX'] <= flux_range[1]):
                    continue

            target = dic['TARGET']
            if targets:
                targets = [x.lower().replace("_", " ") for x in targets]
                target = target.lower().replace("_", " ")
                if target not in targets:
                    continue
            print ("Selected: ", target, date, dic['BAND'], dic['DISP'], dic['DIT'], dic['CATEGORY'], dic['BCDCONFIG'],dic['SEEING'],dic['TAU0'],dic['PWV'],dic['TARGET_FLUX'])
            seeing_value.append(dic['SEEING'])
            tau0_value.append(dic['TAU0'])
            filtered_list_of_dicts.append(dic)
    seeing_value=np.array(seeing_value)
    tau0_value=np.array(tau0_value)
    print("average seeing = {0} +- {1}".format(np.mean(seeing_value),np.std(seeing_value)))
    print("average tau0 = {0} +- {1}".format(np.mean(tau0_value),np.std(tau0_value)))
    return filtered_list_of_dicts


def filter_oi_list_night(
        list_of_dicts: List[Dict],
        dates: Union[str, List[str]] = '2000-01-01',
        bands: Optional[List] =[],
        spectral_resolutions: Optional[List[str]] = [],
        DIT_range: Optional[List[float]] = [],
        targets: Optional[List[str]] = []) -> List[Dict]:
    """Filters the list of oifits-files for one or more nights.

    Parameters
    ----------
    list_of_dicts : list of dict
        List of oifits-files.
    dates : str or list of str
        Date of the observation. The format is ["YYYY-MM-DD"].
        The default is '2000-01-01'.
    bands : list
        List of bands. 'L','M','LM', 'N'. The default is [].
    spectral_resolutions : list of str
        List of spectral resolutions. 'LOW','MED','HIGH'.
        The default is [].
    DIT_range : list of float
        List of DIT-range. The format is [min,max] (s).
        The default is [].
    targets : list of str
        List of targets. The default is [].

    Returns
    -------
    list of dict
        List of filtered oifits-files.
    """
    t=Time(dates)
    mjd0=t.mjd
    filtered_list_of_dicts = []
    if bands:
        bands_new = []
        for i in range(len(bands)):
            if bands[i] == 'M':
                bands_new.append('LM')
            elif bands[i] == 'L':
                bands_new.append('LM')
            else:
                bands_new.append(bands[i])

    for dic in list_of_dicts:
        if dic:
            date = dic['DATEOBS']
            tm=Time(date)
            if (np.abs(tm.mjd-mjd0) > 0.5):
                continue
            if bands:
                if dic['BAND'] not in bands_new:
                    continue
            if spectral_resolutions:
                if dic['DISP'] not in spectral_resolutions:
                    continue
            if DIT_range:
                if not (dic['DIT'] >= DIT_range[0] and dic['DIT'] <= DIT_range[1]):
                    continue

            target = dic['TARGET']
            if targets:
                targets = [x.lower().replace("_", " ") for x in targets]
                target = target.lower().replace("_", " ")
                if target not in targets:
                    continue
            print ("Selected: ", target, date, dic['BAND'], dic['DISP'], dic['DIT'], dic['CATEGORY'])
            filtered_list_of_dicts.append(dic)
    return filtered_list_of_dicts
