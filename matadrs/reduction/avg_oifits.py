import os
from pathlib import Path
from typing import List, Optional
from shutil import copyfile

import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d

import robust


def avg_oifits(infile_list: List[Path], outfile_path: Path,
               headerval: Optional[List[str]] = [],
               avg_baselines: Optional[bool] = False,
               avg_groups: Optional[List[int]] = [6],
               filter_lst_flux: Optional[List] = [],
               filter_lst_corrflux: Optional[List] = [],
               scale_flux: Optional[bool] = False,
               scale_corrflux: Optional[bool] = False,
               scale_wl_range: Optional[List[float]] = [8.3, 9.2],
               avg_func: Optional[str] = 'robustmean') -> None:
    """Averages (.fits)-files.

    Parameters
    ----------
    infile_list : list of pathlib.Path
        List of input files
    outfile_path : Path
        Output file
    headerval : list of str, optional
        List of header keywords
    avg_baselines : bool, optional
        If True, the function will average all the baseline data.
    avg_groups : list of int, optional
        List of group numbers to average.
    filter_lst_flux : Optional[List]
        Filters flux data to average (experimental).
        List of row numbers (starting from 0) to be averaged.
        Example: [[],[0,1],[2,3]]].
    filter_lst_corrflux : list of float, optional
        List of row numbers (starting from 0) to be averaged.
    scale_corrflux : bool, optional
        If True, the function will scale the corr flux data.
    scale_wl_range : list of float, optional
        List of wavelength range to scale.
    avg_func : str, optional
        Function to use for averaging.
        Choices are 'robustmean', 'biweightmean', 'nanmean', 'nanmedian'.
    """
    if avg_func == 'robustmean':
        avgfunc = robust.mean
    if avg_func == 'biweightmean':
        avgfunc = robust.biweightMean
    if avg_func == 'nanmean' or avg_func == 'mean':
        avgfunc = np.nanmean
    if avg_func == 'nanmedian' or avg_func == 'median':
        avgfunc = np.nanmedian

    if os.path.exists(infile_list[0]):
        copyfile(infile_list[0], outfile_path)
    else:
        print(f'ERROR (avg_oifits): File not found: {infile_list[0]}')
        return
    outhdul = fits.open(outfile_path, mode='update')

    inhdul_lst = []
    visamp_lst = []
    visamperr_lst = []
    visphi_lst = []
    visphierr_lst = []
    visamp_ucoord_lst = []
    visamp_vcoord_lst = []
    visamp_sta_index_lst = []
    vis2_lst = []
    vis2err_lst = []
    vis2_ucoord_lst = []
    vis2_vcoord_lst = []
    vis2_sta_index_lst = []
    t3phi_lst = []
    t3phierr_lst = []
    t3phi_u1coord_lst = []
    t3phi_v1coord_lst = []
    t3phi_u2coord_lst = []
    t3phi_v2coord_lst = []
    t3phi_sta_index_lst = []
    flux_lst = []
    fluxerr_lst = []
    wl_lst_flux = []
    wl_lst_vis2 = []
    wl_lst_visamp = []

    j = 0
    for infile in infile_list:
        # NOTE: Read OI_VIS
        if os.path.exists(infile):

            inhdul_lst.append(fits.open(infile, mode='readonly'))
            wl = inhdul_lst[-1]['OI_WAVELENGTH'].data['EFF_WAVE']
            # print(infile)

            # NOTE: Read OI_VIS
            visamp = inhdul_lst[-1]['OI_VIS'].data['VISAMP']
            visamperr = inhdul_lst[-1]['OI_VIS'].data['VISAMPERR']
            visphi = inhdul_lst[-1]['OI_VIS'].data['VISPHI']
            visphierr = inhdul_lst[-1]['OI_VIS'].data['VISPHIERR']
            sta_index = inhdul_lst[-1]['OI_VIS'].data['STA_INDEX']
            ucoord = inhdul_lst[-1]['OI_VIS'].data['UCOORD']
            vcoord = inhdul_lst[-1]['OI_VIS'].data['VCOORD']
            # print('bef', ucoord, sta_index)
            if avg_baselines:
                pbl = np.sqrt(ucoord**2+vcoord**2)
                pbl_sort_idx = np.argsort(pbl)
                jj = 0
                for i in range(len(avg_groups)):
                    sta_idx = sta_index[pbl_sort_idx[jj]]
                    uc = np.nanmean(pbl[pbl_sort_idx[jj:jj+avg_groups[i]]])
                    vc = 0.0
                    for k in range(avg_groups[i]):
                        sta_index[pbl_sort_idx[jj+k]] = sta_idx
                        ucoord[pbl_sort_idx[jj+k]] = uc
                        vcoord[pbl_sort_idx[jj+k]] = vc
                    jj += avg_groups[i]
            # print('aft', ucoord, sta_index)
            for i in range(len(visamp)):
                if np.all(visamp[i] == 0.0):
                    visamp[i] = visamp[i]*np.nan
                    visamperr[i] = visamperr[i]*np.nan
                if np.all(visphi[i] == 0.0):
                    visphi[i] = visphi[i]*np.nan
                    visphierr[i] = visphierr[i]*np.nan
                visamp_lst.append(visamp[i])
                visamperr_lst.append(visamperr[i])
                visphi_lst.append(visphi[i])
                visphierr_lst.append(visphierr[i])
                # if avg_baselines == False:
                visamp_sta_index_lst.append(sta_index[i])
                visamp_ucoord_lst.append(ucoord[i])
                visamp_vcoord_lst.append(vcoord[i])
                # else:
                #     visamp_sta_index_lst.append(sta_index[0])
                #     visamp_ucoord_lst.append(ucoord[0])
                #     visamp_vcoord_lst.append(vcoord[0])
                wl_lst_visamp.append(wl)

            # NOTE: Read OI_VIS2
            vis2 = inhdul_lst[-1]['OI_VIS2'].data['VIS2DATA']
            vis2err = inhdul_lst[-1]['OI_VIS2'].data['VIS2ERR']
            sta_index = inhdul_lst[-1]['OI_VIS2'].data['STA_INDEX']
            ucoord = inhdul_lst[-1]['OI_VIS2'].data['UCOORD']
            vcoord = inhdul_lst[-1]['OI_VIS2'].data['VCOORD']
            if avg_baselines:
                pbl = np.hypot(ucoord, vcoord)
                pbl_sort_idx = np.argsort(pbl)
                jj = 0
                for i in range(len(avg_groups)):
                    sta_idx = sta_index[pbl_sort_idx[jj]]
                    uc = np.nanmean(pbl[pbl_sort_idx[jj:jj+avg_groups[i]]])
                    vc = 0.0
                    for k in range(avg_groups[i]):
                        sta_index[pbl_sort_idx[jj+k]] = sta_idx
                        ucoord[pbl_sort_idx[jj+k]] = uc
                        vcoord[pbl_sort_idx[jj+k]] = vc
                    jj += avg_groups[i]
            for i in range(len(vis2)):
                if np.all(vis2[i] == 0.0):
                    vis2[i] = vis2[i]*np.nan
                    vis2err[i] = vis2err[i]*np.nan
                vis2_lst.append(vis2[i])
                vis2err_lst.append(vis2err[i])
                vis2_sta_index_lst.append(sta_index[i])
                vis2_ucoord_lst.append(ucoord[i])
                vis2_vcoord_lst.append(vcoord[i])
                wl_lst_vis2.append(wl)

            # NOTE: Read OI_T3
            t3phi = inhdul_lst[-1]['OI_T3'].data['T3PHI']
            t3phierr = inhdul_lst[-1]['OI_T3'].data['T3PHIERR']
            sta_index = inhdul_lst[-1]['OI_T3'].data['STA_INDEX']
            u1coord = inhdul_lst[-1]['OI_T3'].data['U1COORD']
            v1coord = inhdul_lst[-1]['OI_T3'].data['V1COORD']
            u2coord = inhdul_lst[-1]['OI_T3'].data['U2COORD']
            v2coord = inhdul_lst[-1]['OI_T3'].data['V2COORD']
            for i in range(len(t3phi)):
                if np.all(t3phi[i] == 0.0):
                    t3phi[i] = t3phi[i]*np.nan
                    t3phierr[i] = t3phierr[i]*np.nan
                t3phi_lst.append(t3phi[i])
                t3phierr_lst.append(t3phierr[i])
                t3phi_sta_index_lst.append(sta_index[i])
                t3phi_u1coord_lst.append(u1coord[i])
                t3phi_v1coord_lst.append(v1coord[i])
                t3phi_u2coord_lst.append(u2coord[i])
                t3phi_v2coord_lst.append(v2coord[i])

            is_flux = True
            # NOTE: Read OI_FLUX
            if len(filter_lst_flux) > 0:
                filter_flux = filter_lst_flux[j]
            else:
                filter_flux = [*range(0, 42, 1)]
            try:
                fluxdata = inhdul_lst[-1]['OI_FLUX'].data['FLUXDATA']
                fluxerr = inhdul_lst[-1]['OI_FLUX'].data['FLUXERR']
                k = 0
                for spectrum, errspectrum in zip(fluxdata, fluxerr):
                    if k in filter_flux:
                        if np.all(spectrum == 0.0):
                            flux_lst.append(spectrum*np.nan)
                            fluxerr_lst.append(errspectrum*np.nan)
                        else:
                            flux_lst.append(spectrum)
                            fluxerr_lst.append(errspectrum)
                        wl_lst_flux.append(wl)
                    k += 1
            except KeyError:
                is_flux = False
                flux_lst.append(np.nan*visamp)
                fluxerr_lst.append(np.nan*visamp)
                wl_lst_flux.append(wl)
        else:
            print(f"WARNING (avg_oifits): File not found: {infile}")
        j += 1

    if not inhdul_lst:
        outhdul.close()
        os.remove(outfile_path)
        print('ERROR (avg_oifits): No files to average.')
        return

    # NOTE: Average fluxes:
    if is_flux:
        flux_arr = np.array(flux_lst)
        fluxerr_arr = np.array(fluxerr_lst)
        # for ii in range(len(flux_arr)):
        #     print(flux_arr[ii])
        # avg_flux = np.nanmean(flux_arr,axis=0)
        avg_flux = avgfunc(flux_arr, axis=0)
        if scale_flux:
            wl_idx = np.logical_and(wl_lst_flux[0]*1e6 > (scale_wl_range[0]),
                                    wl_lst_flux[0]*1e6 <= scale_wl_range[1])
            avg_flux_wl = np.nanmean(avg_flux[wl_idx])
            for i in range(len(flux_lst)):
                flux_lst[i] = flux_lst[i]/np.nanmean(flux_lst[i][wl_idx])*avg_flux_wl
                fluxerr_lst[i] = fluxerr_lst[i]/np.nanmean(flux_lst[i][wl_idx])*avg_flux_wl
            flux_arr = np.array(flux_lst)
            fluxerr_arr = np.array(fluxerr_lst)
            avg_flux = avgfunc(flux_arr, axis=0)

        # avg_flux = np.nanmedian(flux_arr,axis=0)
        # print(avg_flux)
        # print(len(flux_arr))
        if len(flux_arr) > 3:
            # NOTE: Combine two error sources: Standard deviation over the
            # different BCDs, and average error (calculated by the pipeline)
            avg_fluxerr = np.sqrt(np.nanstd(flux_arr, axis=0)**2.0\
                    + np.nanmean(fluxerr_arr, axis=0)**2.0)
        else:
            # WARNING: It may be not the best method for error calculation
            avg_fluxerr = np.nanmean(fluxerr_arr, axis=0)
        outhdul['OI_FLUX'].data = outhdul['OI_FLUX'].data[0:1]
        outhdul['OI_FLUX'].data['FLUXDATA'] = avg_flux
        outhdul['OI_FLUX'].data['FLUXERR'] = avg_fluxerr

    # NOTE: Collect unique station indices from OI_VIS
    sta_index_unique_lst = []
    ucoord_unique_lst = []
    vcoord_unique_lst = []
    sta_index = inhdul_lst[0]['OI_VIS'].data['STA_INDEX']
    sta_index= [list(item) for item in sta_index]
    ucoord = inhdul_lst[0]['OI_VIS'].data['UCOORD']
    vcoord = inhdul_lst[0]['OI_VIS'].data['VCOORD']
    sta_index_unique_lst.append(sta_index[0])
    ucoord_unique_lst.append(ucoord[0])
    vcoord_unique_lst.append(vcoord[0])
    for i in range(1, len(sta_index)):
        if not ((sta_index[i] in sta_index_unique_lst)
                or (sta_index[i][::-1] in sta_index_unique_lst)):
            sta_index_unique_lst.append(sta_index[i])
            ucoord_unique_lst.append(ucoord[i])
            vcoord_unique_lst.append(vcoord[i])

    # NOTE: Average VISAMP and VISPHI
    n_sta_index = len(sta_index_unique_lst)
    outhdul['OI_VIS'].data = outhdul['OI_VIS'].data[0:n_sta_index]
    for k in range(len(sta_index_unique_lst)):
        # NOTE: Collect and average matching visamp data
        visamp_lst_sta = []
        visamperr_lst_sta = []
        visphi_lst_sta = []
        visphierr_lst_sta = []
        for i in range(len(visamp_sta_index_lst)):
            if (((sta_index_unique_lst[k][0] == visamp_sta_index_lst[i][0]) and (sta_index_unique_lst[k][1] == visamp_sta_index_lst[i][1])) \
            or ((sta_index_unique_lst[k][0] == visamp_sta_index_lst[i][1]) and (sta_index_unique_lst[k][1] == visamp_sta_index_lst[i][0]))):
                visamp_lst_sta.append(visamp_lst[i])
                visamperr_lst_sta.append(visamperr_lst[i])
                visphi_lst_sta.append(visphi_lst[i])
                visphierr_lst_sta.append(visphierr_lst[i])

        visamp_arr = np.array(visamp_lst_sta)
        visamperr_arr = np.array(visamperr_lst_sta)
        visphi_arr = np.array(visphi_lst_sta)
        visphierr_arr = np.array(visphierr_lst_sta)
        # print(k)
        # print(visamp_arr)
        avg_visamp = avgfunc(visamp_arr,axis=0)
        if scale_corrflux == True:
            wl_idx = np.logical_and(wl_lst_flux[0]*1e6 > (scale_wl_range[0]),
                                    wl_lst_flux[0]*1e6 <= scale_wl_range[1])
            avg_cflux_wl = np.nanmean(avg_visamp[wl_idx])
            for i in range(len(visamp_lst_sta)):
                visamp_lst_sta[i] = visamp_lst_sta[i]/np.nanmean(visamp_lst_sta[i][wl_idx])*avg_cflux_wl
                visamperr_lst_sta[i] = visamperr_lst_sta[i]/np.nanmean(visamp_lst_sta[i][wl_idx])*avg_cflux_wl
            visamp_arr = np.array(visamp_lst_sta)
            visamperr_arr = np.array(visamperr_lst_sta)
            avg_visamp = avgfunc(visamp_arr, axis=0)
        avg_visphi = np.arctan2(avgfunc(np.sin(visphi_arr*np.pi/180.0), axis=0),
                                avgfunc(np.cos(visphi_arr*np.pi/180.0), axis=0))*180.0/np.pi
        # print(len(visamp_arr))
        if len(visamp_arr) > 3:
            avg_visamperr = np.sqrt(np.nanstd(visamp_arr,axis=0)**2.0 + np.nanmean(visamperr_arr,axis=0)**2.0)
            avg_visphierr = np.sqrt(np.nanstd(visphi_arr,axis=0)**2.0 + np.nanmean(visphierr_arr,axis=0)**2.0)
            # NOTE: Combine two error sources: standard deviation over the different BCDs,
            # and average error (calculated by the pipeline)
        else:
            # WARNING: It may be not the best method for error calculation
            avg_visamperr = np.nanmean(visamperr_arr, axis=0)
            avg_visphierr = np.nanmean(visphierr_arr, axis=0)
        # print(avg_visamp)
        # print(outhdul['OI_VIS'].data['VISAMP'][k])
        outhdul['OI_VIS'].data['VISAMP'][k] = avg_visamp
        outhdul['OI_VIS'].data['VISAMPERR'][k] = avg_visamperr
        outhdul['OI_VIS'].data['VISPHI'][k] = avg_visphi
        outhdul['OI_VIS'].data['VISPHIERR'][k] = avg_visphierr
        outhdul['OI_VIS'].data['STA_INDEX'][k] = sta_index_unique_lst[k]
        outhdul['OI_VIS'].data['UCOORD'][k] = ucoord_unique_lst[k]
        outhdul['OI_VIS'].data['VCOORD'][k] = vcoord_unique_lst[k]

    # NOTE: Collect unique station indices from OI_VIS2
    sta_index_unique_lst = []
    ucoord_unique_lst = []
    vcoord_unique_lst = []
    sta_index = inhdul_lst[0]['OI_VIS2'].data['STA_INDEX']
    sta_index = [list(item) for item in sta_index]
    ucoord = inhdul_lst[0]['OI_VIS2'].data['UCOORD']
    vcoord = inhdul_lst[0]['OI_VIS2'].data['VCOORD']
    sta_index_unique_lst.append(sta_index[0])
    ucoord_unique_lst.append(ucoord[0])
    vcoord_unique_lst.append(vcoord[0])
    for i in range(1, len(sta_index)):
        if not ((sta_index[i] in sta_index_unique_lst)
                or (sta_index[i][::-1] in sta_index_unique_lst)):
            sta_index_unique_lst.append(sta_index[i])
            ucoord_unique_lst.append(ucoord[i])
            vcoord_unique_lst.append(vcoord[i])

    # NOTE: Average VIS2
    n_sta_index = len(sta_index_unique_lst)
    outhdul['OI_VIS2'].data = outhdul['OI_VIS2'].data[0:n_sta_index]
    for k in range(len(sta_index_unique_lst)):
        # NOTE: Collect and average matching vis2 data
        vis2_lst_sta = []
        vis2err_lst_sta = []
        for i in range(len(vis2_sta_index_lst)):
            if (((sta_index_unique_lst[k][0] == vis2_sta_index_lst[i][0]) and (sta_index_unique_lst[k][1] == vis2_sta_index_lst[i][1])) \
            or ((sta_index_unique_lst[k][0] == vis2_sta_index_lst[i][1]) and (sta_index_unique_lst[k][1] == vis2_sta_index_lst[i][0]))):
                vis2_lst_sta.append(vis2_lst[i])
                vis2err_lst_sta.append(vis2err_lst[i])
        vis2_arr = np.array(vis2_lst_sta)
        vis2err_arr = np.array(vis2err_lst_sta)
        avg_vis2 = avgfunc(vis2_arr, axis=0)
        if len(vis2_arr) > 3:
            # NOTE: Combine two error sources: standard deviation over the different BCDs,
            # and average error (calculated by the pipeline)
            avg_vis2err = np.sqrt(np.nanstd(vis2_arr, axis=0)**2.0\
                    + np.nanmean(vis2err_arr, axis=0)**2.0)
        else:
            # WARNING: It may be not the best method for error calculation
            avg_vis2err = np.nanmean(vis2err_arr, axis=0)
        outhdul['OI_VIS2'].data['VIS2DATA'][k] = avg_vis2
        outhdul['OI_VIS2'].data['VIS2ERR'][k] = avg_vis2err
        outhdul['OI_VIS2'].data['STA_INDEX'][k] = sta_index_unique_lst[k]
        outhdul['OI_VIS2'].data['UCOORD'][k] = ucoord_unique_lst[k]
        outhdul['OI_VIS2'].data['VCOORD'][k] = vcoord_unique_lst[k]

    # NOTE: Collect unique station indices from OI_T3
    sta_index_unique_lst = []
    sta_index_unique_lst_sorted = []
    u1coord_unique_lst = []
    v1coord_unique_lst = []
    u2coord_unique_lst = []
    v2coord_unique_lst = []
    sta_index = inhdul_lst[0]['OI_T3'].data['STA_INDEX']
    sta_index = [list(item) for item in sta_index]
    u1coord = inhdul_lst[0]['OI_T3'].data['U1COORD']
    v1coord = inhdul_lst[0]['OI_T3'].data['V1COORD']
    u2coord = inhdul_lst[0]['OI_T3'].data['U2COORD']
    v2coord = inhdul_lst[0]['OI_T3'].data['V2COORD']
    sta_index_unique_lst.append(sta_index[0])
    sta_index_unique_lst_sorted.append(sorted(sta_index[0]))
    # print(sorted(sta_index[0]))
    u1coord_unique_lst.append(u1coord[0])
    v1coord_unique_lst.append(v1coord[0])
    u2coord_unique_lst.append(u2coord[0])
    v2coord_unique_lst.append(v2coord[0])
    for i in range(1, len(sta_index)):
        if not ((sorted(sta_index[i]) in sta_index_unique_lst_sorted)):
            sta_index_unique_lst.append(sta_index[i])
            sta_index_unique_lst_sorted.append(sorted(sta_index[i]))
            u1coord_unique_lst.append(u1coord[i])
            v1coord_unique_lst.append(v1coord[i])
            u2coord_unique_lst.append(u2coord[i])
            v2coord_unique_lst.append(v2coord[i])

    # NOTE: Average T3PHI
    n_sta_index = len(sta_index_unique_lst)
    outhdul['OI_T3'].data = outhdul['OI_T3'].data[0:n_sta_index]
    for k in range(len(sta_index_unique_lst)):
        # NOTE: Collect and average matching vis2 data
        t3phi_lst_sta = []
        t3phierr_lst_sta = []
        # print('k',k,sta_index_unique_lst_sorted[k])
        for i in range(len(t3phi_sta_index_lst)):
            # print('i',i,sorted(t3phi_sta_index_lst[i]))
            if sta_index_unique_lst_sorted[k] == sorted(t3phi_sta_index_lst[i]):
                t3phi_lst_sta.append(t3phi_lst[i])
                t3phierr_lst_sta.append(t3phierr_lst[i])
        t3phi_arr = np.array(t3phi_lst_sta)
        t3phierr_arr = np.array(t3phierr_lst_sta)
        # avg_t3phi = np.nanmean(t3phi_arr,axis=0)
        avg_t3phi = np.arctan2(avgfunc(np.sin(t3phi_arr*np.pi/180.0),axis=0),avgfunc(np.cos(t3phi_arr*np.pi/180.0),axis=0))*180.0/np.pi
        if len(t3phi_arr) > 3:
            # NOTE: Combine two error sources: standard deviation over the different BCDs,
            # and average error (calculated by the pipeline)
            avg_t3phierr = np.sqrt(np.nanstd(t3phi_arr,axis=0)**2.0 + np.nanmean(t3phierr_arr,axis=0)**2.0)
        else:
            # WARNING: It may be not the best method for error calculation
            avg_t3phierr = np.nanmean(t3phierr_arr, axis=0)
        outhdul['OI_T3'].data['T3PHI'][k] = avg_t3phi
        outhdul['OI_T3'].data['T3PHIERR'][k] = avg_t3phierr
        outhdul['OI_T3'].data['STA_INDEX'][k] = sta_index_unique_lst[k]
        outhdul['OI_T3'].data['U1COORD'][k] = u1coord_unique_lst[k]
        outhdul['OI_T3'].data['V1COORD'][k] = v1coord_unique_lst[k]
        outhdul['OI_T3'].data['U2COORD'][k] = u2coord_unique_lst[k]
        outhdul['OI_T3'].data['V2COORD'][k] = v2coord_unique_lst[k]

    for dic in headerval:
        outhdul[0].header[dic['key']] = dic['value']

    # NOTE: Changes are written back to original fits
    outhdul.flush()
    outhdul.close()
    for inhdul in inhdul_lst:
        inhdul.close()


# oi_types_list = [ ['vis2','t3'], ['visamp'] ]
def oifits_patchwork(
        infile_list: List[Path], outfile_path: Path,
        oi_types_list: Optional[List[List[str]]] = [['vis2', 'visamp', 'visphi', 't3', 'flux']],
        headerval: Optional[List[str]] = []) -> None:
    """

    Parameters
    ----------
    infile_list : list of pathlib.Path
        List of input files.
    outfile_path : pathlib.Path
        Output file.
    oi_types_list : list of list of str
        List of OI types.
    headerval : list of str
        List of header keys.
    """
    # print(infile_list)
    if os.path.exists(infile_list[0]):
        copyfile(infile_list[0], outfile_path)
    else:
        print(f'ERROR (oifits_patchwork): File not found: {infile_list[0]}')
        return
    outhdul = fits.open(outfile_path, mode='update')

    n_oi_types_list = len(oi_types_list)
    for i in range(n_oi_types_list):
        # print(i)
        oi_types = oi_types_list[i]
        infile = infile_list[i]
        # print(infile)
        inhdul = fits.open(infile, mode='readonly')

        for oi_type in oi_types:
            if oi_type == 'vis2':
                outhdul['OI_VIS2'].data = inhdul['OI_VIS2'].data
            if oi_type == 't3':
                # n_data = len(inhdul['OI_T3'].data)
                # print('T3Phi',n_data)
                # outhdul['OI_T3'].data = outhdul['OI_T3'].data[0:n_data]
                outhdul['OI_T3'].data = inhdul['OI_T3'].data
            if oi_type == 'visamp':
                try:
                    outhdul[0].header['HIERARCH ESO PRO CAL NAME'] = inhdul[0].header['HIERARCH ESO PRO CAL NAME']
                    outhdul[0].header['HIERARCH ESO PRO CAL RA'] = inhdul[0].header['HIERARCH ESO PRO CAL RA']  
                    outhdul[0].header['HIERARCH ESO PRO CAL DEC'] = inhdul[0].header['HIERARCH ESO PRO CAL DEC'] 
                    outhdul[0].header['HIERARCH ESO PRO CAL AIRM'] = inhdul[0].header['HIERARCH ESO PRO CAL AIRM']
                    outhdul[0].header['HIERARCH ESO PRO CAL FWHM'] = inhdul[0].header['HIERARCH ESO PRO CAL FWHM']
                    outhdul[0].header['HIERARCH ESO PRO CAL TAU0'] = inhdul[0].header['HIERARCH ESO PRO CAL TAU0']
                    outhdul[0].header['HIERARCH ESO PRO CAL TPL START'] = inhdul[0].header['HIERARCH ESO PRO CAL TPL START']     
                    outhdul[0].header['HIERARCH ESO PRO CAL DB NAME'] = inhdul[0].header['HIERARCH ESO PRO CAL DB NAME']    
                    outhdul[0].header['HIERARCH ESO PRO CAL DB DBNAME'] = inhdul[0].header['HIERARCH ESO PRO CAL DB DBNAME']  
                    outhdul[0].header['HIERARCH ESO PRO CAL DB RA'] = inhdul[0].header['HIERARCH ESO PRO CAL DB RA']      
                    outhdul[0].header['HIERARCH ESO PRO CAL DB DEC'] = inhdul[0].header['HIERARCH ESO PRO CAL DB DEC']     
                    outhdul[0].header['HIERARCH ESO PRO CAL DB DIAM'] = inhdul[0].header['HIERARCH ESO PRO CAL DB DIAM']    
                    outhdul[0].header['HIERARCH ESO PRO CAL DB ERRDIAM'] = inhdul[0].header['HIERARCH ESO PRO CAL DB ERRDIAM'] 
                    outhdul[0].header['HIERARCH ESO PRO CAL DB SEPARATION'] = inhdul[0].header['HIERARCH ESO PRO CAL DB SEPARATION'] 
                except KeyError as e:
                    print(e)

                outhdul['OI_VIS'].header['AMPTYP'] = inhdul['OI_VIS'].header['AMPTYP']
                outhdul['OI_VIS'].data = inhdul['OI_VIS'].data
                # NOTE: Look up visphi
                for j in range(n_oi_types_list):
                    if 'visphi' in oi_types_list[j]:
                        infile2 = infile_list[j]
                        inhdul2 = fits.open(infile2, mode='readonly')
                        visphi = inhdul2['OI_VIS'].data['VISPHI']
                        visphierr = inhdul2['OI_VIS'].data['VISPHIERR']

                        # NOTE: Match station indices
                        sta_index_visamp = outhdul['OI_VIS'].data['STA_INDEX']
                        sta_index_visamp = [ list(item) for item in sta_index_visamp ]
                        sta_index_visphi = inhdul2['OI_VIS'].data['STA_INDEX']
                        sta_index_visphi = [ list(item) for item in sta_index_visphi ]
                        for k in range(len(sta_index_visamp)):
                            for l in range(len(sta_index_visphi)):
                                if ( (sta_index_visamp[k] == sta_index_visphi[l]) \
                                    or (sta_index_visamp[k][::-1] == sta_index_visphi[l] )):
                                    outhdul['OI_VIS'].data['VISPHI'][k] = inhdul2['OI_VIS'].data['VISPHI'][l]
                                    outhdul['OI_VIS'].data['VISPHIERR'][k] = inhdul2['OI_VIS'].data['VISPHIERR'][l]
            if oi_type == 'flux':
                try:
                    outhdul['OI_FLUX'].data = inhdul['OI_FLUX'].data
                except KeyError as e:
                    pass

    for dic in headerval:
        del outhdul[0].header[dic['key']]
        outhdul[0].header[dic['key']] = dic['value']

    # NOTE: Changes are written back to original.fits
    outhdul.flush()
    outhdul.close()
    inhdul.close()
    inhdul2.close()


def calc_vis_from_corrflux(input_corrflux_file,input_totalflux_file,outfile_path):
    copyfile(input_corrflux_file, outfile_path)
    outhdul  = fits.open(outfile_path, mode='update')

    inhdul_corr = fits.open(input_corrflux_file, mode='readonly')
    inhdul_tot = fits.open(input_totalflux_file, mode='readonly')

    # NOTE: Read total spectra
    wl_flux = inhdul_tot['OI_WAVELENGTH'].data['EFF_WAVE']
    flux = inhdul_tot['OI_FLUX'].data['FLUXDATA'][0]
    fluxerr = inhdul_tot['OI_FLUX'].data['FLUXERR'][0]
    # print(flux,fluxerr)

    # NOTE: Read correlated spectra
    corrflux = inhdul_corr['OI_VIS'].data['VISAMP']
    corrfluxerr = inhdul_corr['OI_VIS'].data['VISAMPERR']
    wl_corrflux = inhdul_corr['OI_WAVELENGTH'].data['EFF_WAVE']
    # print(corrflux,corrfluxerr)

    if not len(outhdul['OI_VIS'].data['VISAMP'][0]) == len(flux):
        # NOTE: Interpolate the flux data to the wavelengths of the correlated flux
        f = interp1d(wl_flux, flux, kind='cubic')
        flux_resamp = f(wl_corrflux)
        f = interp1d(wl_flux, fluxerr, kind='cubic')
        fluxerr_resamp = f(wl_corrflux)
        flux = flux_resamp
        fluxerr = fluxerr_resamp
        outhdul['OI_FLUX'].data['FLUXDATA'] = flux
        outhdul['OI_FLUX'].data['FLUXERR'] = fluxerr
        
    for k in range(len(outhdul['OI_VIS'].data['VISAMP'])):
        # NOTE: Collect and average matching vis2 data
        vis = (corrflux[k]/flux)
        viserr = vis*np.sqrt((corrfluxerr[k]/corrflux[k])**2 + (fluxerr/flux)**2)
        vis2 = vis**2.0
        vis2err = 2.0*vis*viserr
        # print(viserr,vis2err)

        # WARNING: this is not how errors should be calculated
        outhdul['OI_VIS2'].data['VIS2DATA'][k] = vis2
        outhdul['OI_VIS2'].data['VIS2ERR'][k] = vis2err
        outhdul['OI_VIS2'].data['STA_INDEX'][k] = inhdul_corr['OI_VIS'].data['STA_INDEX'][k]
        outhdul['OI_VIS2'].data['UCOORD'][k] = inhdul_corr['OI_VIS'].data['UCOORD'][k]
        outhdul['OI_VIS2'].data['VCOORD'][k] = inhdul_corr['OI_VIS'].data['VCOORD'][k]

    # NOTE: Changes are written back to original.fits
    outhdul.flush()
    outhdul.close()
    inhdul_corr.close()
    inhdul_tot.close()


def oifits_restrict_wavelengths(
        infile_path: Path, outfile_path: Path,
        sel_wl: Optional[float] = 3.0,
        bandwidth: Optional[float] = 1000.0,
        wl_range: Optional[List[float]] = [],
        do_flag: Optional[bool] = False) -> None:
    """Restricts the wavelengths of an OIFITS file.

    Parameters
    ----------
    infile_path : pathlib.Path
        Input file.
    outfile_path : pathlib.Path
        Output file.
    sel_wl : float, optional
        The selected wavelength.
    bandwidth : float, optional
        The bandwidth of the desired wavelength.
    wl_range : list of float, optional
        The wavelength range.
    do_flag : bool, optional
        If True, flag the wavelengths we do not want to use.
    """
    copyfile(infile_path, outfile_path)
    outhdul = fits.open(outfile_path, mode='update')

    hdu = fits.open(infile_path, mode='readonly')
    try:
        wl = hdu['OI_WAVELENGTH'].data['EFF_WAVE']
        band = hdu['OI_WAVELENGTH'].data['EFF_BAND']
    except KeyError as e:
        print("ERROR: No OI_WAVELENGTH table! ")
        print(e)
        return

    if wl_range:
        idx = np.logical_and(wl*1e6 > (wl_range[0]), wl*1e6 <= wl_range[1])
    else:
        idx = np.logical_or(wl*1e6 < (sel_wl-bandwidth/2.0),wl*1e6 >= (sel_wl+bandwidth/2.0))
    # wl_new = wl[idx]
    # band_new = band[idx]

    if do_flag:
        # NOTE: Flag the wavelengths we do not want to use
        try:
            # hdu['OI_VIS'].data['VISAMP']
            # hdu['OI_VIS'].data['VISAMPERR']
            # hdu['OI_VIS'].data['VISPHI']
            # hdu['OI_VIS'].data['VISPHIERR']
            # hdu['OI_VIS'].data['FLAG']
            for k in range(len(hdu['OI_VIS'].data['VISAMP'])):
                outhdul['OI_VIS'].data['FLAG'][k][idx] = True
        except KeyError as e:
            print("WARNING: No OI_VIS table!")
            print(e)

        try:
            for k in range(len(hdu['OI_VIS2'].data['VIS2DATA'])):
                # outhdul['OI_VIS2'].data['VIS2DATA'][k][idx] = np.nan
                # outhdul['OI_VIS2'].data['VIS2ERR'][k][idx] = np.nan
                outhdul['OI_VIS2'].data['FLAG'][k][idx] = True
        except KeyError as e:
            print("WARNING: No OI_VIS2 table!")
            print(e)

        try:
            # hdu['TF2'].data['TF2']
            # hdu['TF2'].data['TF2ERR']
            for k in range(len(hdu['TF2'].data['TF2'])):
                outhdul['TF2'].data['FLAG'][k][idx] = True
        except KeyError as e:
            pass
            # print("WARNING: No OI_TF2 table!")
            # print(e)

        try:
            # hdu['OI_T3'].data['T3AMP']
            # hdu['OI_T3'].data['T3AMPERR']
            # hdu['OI_T3'].data['T3PHI']
            # hdu['OI_T3'].data['T3PHIERR']
            for k in range(len(hdu['OI_T3'].data['T3AMP'])):
                outhdul['OI_T3'].data['FLAG'][k][idx] = True
        except KeyError as e:
            print("WARNING: No OI_T3 table!")
            print(e)

        try:
            # hdu['OI_FLUX'].data['FLUXDATA']
            # hdu['OI_FLUX'].data['FLUXERR']
            for k in range(len(hdu['OI_FLUX'].data['FLUXDATA'])):
                outhdul['OI_FLUX'].data['FLAG'][k][idx] = True
        except KeyError as e:
            pass
            # print("WARNING: No OI_FLUX table!")
            # print(e)
    else:
        new_wave = hdu['OI_WAVELENGTH'].data['EFF_WAVE'][idx]
        outhdul['OI_WAVELENGTH'].data = np.resize(outhdul['OI_WAVELENGTH'].data,(len(new_wave),))
        outhdul['OI_WAVELENGTH'].data['EFF_WAVE'] = new_wave
        outhdul['OI_WAVELENGTH'].data['EFF_BAND'] = hdu['OI_WAVELENGTH'].data['EFF_BAND'][idx]

        Nwl = len(new_wave)
        if 'OI_VIS2' in hdu:
            n_rows = len(outhdul['OI_VIS2'].data['VIS2DATA'])
            # OI_VIS2
            NV2 = len(outhdul['OI_VIS2'].data['VIS2DATA'])
            del outhdul['OI_VIS2']
            c1 = fits.Column(name='TARGET_ID', format='I')
            c2 = fits.Column(name='TIME', format='D', unit='s')
            c3 = fits.Column(name='MJD', format='D', unit='day')
            c4 = fits.Column(name='INT_TIME', format='D', unit='s')
            c5 = fits.Column(name='VIS2DATA', format='%dD' % Nwl)
            c6 = fits.Column(name='VIS2ERR', format='%dD' % Nwl)
            c7 = fits.Column(name='UCOORD', format='D', unit='m')
            c8 = fits.Column(name='VCOORD', format='D', unit='m')
            c9 = fits.Column(name='STA_INDEX', format='2I')
            c10 = fits.Column(name='FLAG', format='%dL' % Nwl)
            oi_vis2_hdu = fits.BinTableHDU.from_columns(
                    [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10],
                    name='OI_VIS2', nrows=NV2, fill=True)
            outhdul.append(oi_vis2_hdu)

            outhdul['OI_VIS2'].header['DATE-OBS'] = hdu['OI_VIS2'].header['DATE-OBS']
            outhdul['OI_VIS2'].header['EXTNAME'] = hdu['OI_VIS2'].header['EXTNAME']
            outhdul['OI_VIS2'].header['EXTVER'] = hdu['OI_VIS2'].header['EXTVER']
            outhdul['OI_VIS2'].header['OI_REVN'] = hdu['OI_VIS2'].header['OI_REVN']
            outhdul['OI_VIS2'].header['ARRNAME'] = hdu['OI_VIS2'].header['ARRNAME']
            outhdul['OI_VIS2'].header['INSNAME'] = hdu['OI_VIS2'].header['INSNAME']
            outhdul['OI_VIS2'].data['TARGET_ID'] = hdu['OI_VIS2'].data['TARGET_ID']
            outhdul['OI_VIS2'].data['TIME'] = hdu['OI_VIS2'].data['TIME']
            outhdul['OI_VIS2'].data['MJD'] = hdu['OI_VIS2'].data['MJD']
            outhdul['OI_VIS2'].data['INT_TIME'] = hdu['OI_VIS2'].data['INT_TIME']
            outhdul['OI_VIS2'].data['UCOORD'] = hdu['OI_VIS2'].data['UCOORD']
            outhdul['OI_VIS2'].data['VCOORD'] = hdu['OI_VIS2'].data['VCOORD']
            outhdul['OI_VIS2'].data['STA_INDEX'] = hdu['OI_VIS2'].data['STA_INDEX']
            for j in range(NV2):
                outhdul['OI_VIS2'].data['VIS2DATA'][j] = hdu['OI_VIS2'].data['VIS2DATA'][j][idx]
                outhdul['OI_VIS2'].data['FLAG'][j] = hdu['OI_VIS2'].data['FLAG'][j][idx]
                outhdul['OI_VIS2'].data['VIS2ERR'][j] = hdu['OI_VIS2'].data['VIS2ERR'][j][idx]

        if 'OI_T3' in hdu:
            # OI_T3
            NT3 = len(hdu['OI_T3'].data['T3PHI'])
            del outhdul['OI_T3']
            c1 = fits.Column(name='TARGET_ID', format='1I')
            c2 = fits.Column(name='TIME', format='1D', unit='s')
            c3 = fits.Column(name='MJD', format='1D', unit='day')
            c4 = fits.Column(name='INT_TIME', format='1D', unit='s')
            c5 = fits.Column(name='T3AMP', format='%dD' % Nwl)
            c6 = fits.Column(name='T3AMPERR', format='%dD' % Nwl)
            c7 = fits.Column(name='T3PHI', format='%dD' % Nwl, unit='deg')
            c8 = fits.Column(name='T3PHIERR', format='%dD' % Nwl, unit='deg')
            c9 = fits.Column(name='U1COORD', format='1D', unit='m')
            c10 = fits.Column(name='V1COORD', format='1D', unit='m')
            c11 = fits.Column(name='U2COORD', format='1D', unit='m')
            c12 = fits.Column(name='V2COORD', format='1D', unit='m')
            c13 = fits.Column(name='STA_INDEX', format='3I')
            c14 = fits.Column(name='FLAG', format='%dL' % Nwl)
            oi_t3_hdu = fits.BinTableHDU.from_columns(
                    [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14],
                    name='OI_T3', nrows=NT3, fill=True)
            outhdul.append(oi_t3_hdu)

            outhdul['OI_T3'].header['DATE-OBS']= hdu['OI_T3'].header['DATE-OBS']
            outhdul['OI_T3'].header['EXTNAME'] = hdu['OI_T3'].header['EXTNAME']
            outhdul['OI_T3'].header['EXTVER'] = hdu['OI_T3'].header['EXTVER']
            outhdul['OI_T3'].header['OI_REVN'] = hdu['OI_T3'].header['OI_REVN']
            outhdul['OI_T3'].header['ARRNAME'] = hdu['OI_T3'].header['ARRNAME']
            outhdul['OI_T3'].header['INSNAME'] = hdu['OI_T3'].header['INSNAME']
            outhdul['OI_T3'].data['TARGET_ID'] = hdu['OI_T3'].data['TARGET_ID']
            outhdul['OI_T3'].data['TIME'] = hdu['OI_T3'].data['TIME']
            outhdul['OI_T3'].data['MJD'] = hdu['OI_T3'].data['MJD']
            outhdul['OI_T3'].data['INT_TIME'] = hdu['OI_T3'].data['INT_TIME']
            outhdul['OI_T3'].data['U1COORD'] = hdu['OI_T3'].data['U1COORD']
            outhdul['OI_T3'].data['V1COORD'] = hdu['OI_T3'].data['V1COORD']
            outhdul['OI_T3'].data['U2COORD'] = hdu['OI_T3'].data['U2COORD']
            outhdul['OI_T3'].data['V2COORD'] = hdu['OI_T3'].data['V2COORD']
            outhdul['OI_T3'].data['STA_INDEX'] = hdu['OI_T3'].data['STA_INDEX']
            for j in range(NT3):
                outhdul['OI_T3'].data['T3AMP'][j] = hdu['OI_T3'].data['T3AMP'][j][idx]
                outhdul['OI_T3'].data['T3AMPERR'][j] = hdu['OI_T3'].data['T3AMPERR'][j][idx]
                outhdul['OI_T3'].data['T3PHI'][j] = hdu['OI_T3'].data['T3PHI'][j][idx]
                outhdul['OI_T3'].data['T3PHIERR'][j] = hdu['OI_T3'].data['T3PHIERR'][j][idx]
                outhdul['OI_T3'].data['FLAG'][j] = hdu['OI_T3'].data['FLAG'][j][idx]

        if 'OI_VIS' in hdu:
            # OI_VIS
            NV = len(outhdul['OI_VIS'].data['VISAMP'])
            del outhdul['OI_VIS']
            c1 = fits. Column(name='TARGET_ID', format='1I')
            c2 = fits. Column(name='TIME', format='1D', unit='s')
            c3 = fits. Column(name='MJD', format='1D', unit='day')
            c4 = fits. Column(name='INT_TIME', format='1D',unit='s')
            c5 = fits. Column(name='VISAMP', format='%dD' % Nwl)
            c6 = fits. Column(name='VISAMPERR', format='%dD' % Nwl)
            c7 = fits. Column(name='VISPHI', format='%dD' % Nwl, unit='deg')
            c8 = fits. Column(name='VISPHIERR', format='%dD' % Nwl, unit='deg')
            c9 = fits. Column(name='UCOORD', format='1D', unit='m')
            c10= fits. Column(name='VCOORD', format='1D', unit='m')
            c11= fits. Column(name='STA_INDEX', format='2I')
            c12= fits. Column(name='FLAG', format='%dL' % Nwl)
            oi_vis_hdu = fits.BinTableHDU.from_columns(
                    [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12],
                    name='OI_VIS', nrows=NV, fill=True)
            outhdul.append(oi_vis_hdu)

            outhdul['OI_VIS'].header['EXTNAME'] = hdu['OI_VIS'].header['EXTNAME']
            outhdul['OI_VIS'].header['EXTVER'] = hdu['OI_VIS'].header['EXTVER']
            outhdul['OI_VIS'].header['OI_REVN'] = hdu['OI_VIS'].header['OI_REVN']
            outhdul['OI_VIS'].header['ARRNAME'] = hdu['OI_VIS'].header['ARRNAME']
            outhdul['OI_VIS'].header['INSNAME'] = hdu['OI_VIS'].header['INSNAME']
            outhdul['OI_VIS'].header['DATE-OBS']= hdu['OI_VIS'].header['DATE-OBS']
            outhdul['OI_VIS'].header['AMPTYP'] = hdu['OI_VIS'].header['AMPTYP']
            outhdul['OI_VIS'].header['PHITYP'] = hdu['OI_VIS'].header['PHITYP']
            outhdul['OI_VIS'].data['TARGET_ID'] = hdu['OI_VIS'].data['TARGET_ID']
            outhdul['OI_VIS'].data['TIME'] = hdu['OI_VIS'].data['TIME']
            outhdul['OI_VIS'].data['MJD'] = hdu['OI_VIS'].data['MJD']
            outhdul['OI_VIS'].data['INT_TIME'] = hdu['OI_VIS'].data['INT_TIME']
            outhdul['OI_VIS'].data['UCOORD'] = hdu['OI_VIS'].data['UCOORD']
            outhdul['OI_VIS'].data['VCOORD'] = hdu['OI_VIS'].data['VCOORD']
            outhdul['OI_VIS'].data['STA_INDEX'] = hdu['OI_VIS'].data['STA_INDEX']
            for j in range(NV):
                outhdul['OI_VIS'].data['VISAMP'][j] = hdu['OI_VIS'].data['VISAMP'][j][idx]
                outhdul['OI_VIS'].data['VISAMPERR'][j] = hdu['OI_VIS'].data['VISAMPERR'][j][idx]
                outhdul['OI_VIS'].data['VISPHI'][j] = hdu['OI_VIS'].data['VISPHI'][j][idx]
                outhdul['OI_VIS'].data['VISPHIERR'][j] = hdu['OI_VIS'].data['VISPHIERR'][j][idx]
                outhdul['OI_VIS'].data['FLAG'][j] = hdu['OI_VIS'].data['FLAG'][j][idx]

        if 'TF2' in hdu:
            NV2 = len(hdu['TF2'].data['TF2'])
            del outhdul['TF2']
            c1 = fits. Column(name='TIME', format='1D', unit='s')
            c2 = fits. Column(name='MJD', format='1D', unit='day')
            c3 = fits. Column(name='INT_TIME', format='1D', unit='s')
            c4 = fits. Column(name='TF2', format='%dD' % Nwl)
            c5 = fits. Column(name='TF2ERR', format='%dD' % Nwl)
            c6 = fits. Column(name='TF', format='%dD' % Nwl)
            c7 = fits. Column(name='TFERR', format='%dD' % Nwl)
            c8 = fits. Column(name='STA_INDEX', format='2I')
            oi_tf2_hdu = fits.BinTableHDU.from_columns(
                    [c1, c2, c3, c4, c5, c6, c7, c8], name='TF2', nrows=NV2, fill=True)
            outhdul.append(oi_tf2_hdu)

            outhdul['TF2'].header['DATE-OBS']= hdu['TF2'].header['DATE-OBS']
            outhdul['TF2'].header['EXTNAME'] = hdu['TF2'].header['EXTNAME']
            outhdul['TF2'].header['EXTVER']  = hdu['TF2'].header['EXTVER']
            outhdul['TF2'].header['OI_REVN'] = hdu['TF2'].header['OI_REVN']
            outhdul['TF2'].header['ARRNAME'] = hdu['TF2'].header['ARRNAME']
            outhdul['TF2'].header['INSNAME'] = hdu['TF2'].header['INSNAME']
            outhdul['TF2'].data['TIME'] = hdu['TF2'].data['TIME']
            outhdul['TF2'].data['MJD'] = hdu['TF2'].data['MJD']
            outhdul['TF2'].data['INT_TIME'] = hdu['TF2'].data['INT_TIME']
            outhdul['TF2'].data['STA_INDEX'] = hdu['TF2'].data['STA_INDEX']
            for j in range(NV2):
                outhdul['TF2'].data['TF2'][j] = hdu['TF2'].data['TF2'][j][idx]
                outhdul['TF2'].data['TF'][j] = hdu['TF2'].data['TF'][j][idx]
                outhdul['TF2'].data['TF2ERR'][j]= hdu['TF2'].data['TF2ERR'][j][idx]
                outhdul['TF2'].data['TFERR'][j] = hdu['TF2'].data['TFERR'][j][idx]

        if 'OI_FLUX' in hdu:
            # OI_FLUX
            NF = len(hdu['OI_FLUX'].data['FLUXDATA'])
            del outhdul['OI_FLUX']
            c1 = fits. Column(name='TARGET_ID', format='1I')
            c2 = fits. Column(name='TIME', format='1D', unit='s')
            c3 = fits. Column(name='MJD', format='1D', unit='day')
            c4 = fits. Column(name='INT_TIME', format='1D', unit='s')
            c5 = fits. Column(name='FLUXDATA', format='%dD' % Nwl)
            c6 = fits. Column(name='FLUXERR', format='%dD' % Nwl)
            c7 = fits. Column(name='STA_INDEX', format='1I')
            c8 = fits. Column(name='FLAG', format='%dL' % Nwl)
            oi_flux_hdu = fits.BinTableHDU.from_columns(
                    [c1, c2, c3, c4, c5, c6, c7, c8], name='OI_FLUX', nrows=NF, fill=True)
            outhdul.append(oi_flux_hdu)

            outhdul['OI_FLUX'].header['EXTNAME'] = hdu['OI_FLUX'].header['EXTNAME']
            outhdul['OI_FLUX'].header['EXTVER'] = hdu['OI_FLUX'].header['EXTVER']
            outhdul['OI_FLUX'].header['OI_REVN'] = hdu['OI_FLUX'].header['OI_REVN']
            outhdul['OI_FLUX'].header['DATE-OBS']= hdu['OI_FLUX'].header['DATE-OBS']
            outhdul['OI_FLUX'].header['ARRNAME'] = hdu['OI_FLUX'].header['ARRNAME']
            outhdul['OI_FLUX'].header['INSNAME'] = hdu['OI_FLUX'].header['INSNAME']
            outhdul['OI_FLUX'].header['FOV'] = hdu['OI_FLUX'].header['FOV']
            outhdul['OI_FLUX'].header['FOVTYPE'] = hdu['OI_FLUX'].header['FOVTYPE']
            outhdul['OI_FLUX'].header['CALSTAT'] = hdu['OI_FLUX'].header['CALSTAT']
            outhdul['OI_FLUX'].data['TARGET_ID'] = hdu['OI_FLUX'].data['TARGET_ID']
            outhdul['OI_FLUX'].data['TIME'] = hdu['OI_FLUX'].data['TIME']
            outhdul['OI_FLUX'].data['MJD'] = hdu['OI_FLUX'].data['MJD']
            outhdul['OI_FLUX'].data['INT_TIME'] = hdu['OI_FLUX'].data['INT_TIME']
            outhdul['OI_FLUX'].data['STA_INDEX'] = hdu['OI_FLUX'].data['STA_INDEX']
            for j in range(NF):
                outhdul['OI_FLUX'].data['FLUXDATA'][j] = hdu['OI_FLUX'].data['FLUXDATA'][j][idx]
                outhdul['OI_FLUX'].data['FLUXERR'][j] = hdu['OI_FLUX'].data['FLUXERR'][j][idx]
                outhdul['OI_FLUX'].data['FLAG'][j] = hdu['OI_FLUX'].data['FLAG'][j][idx]

    # NOTE: Changes are written back to original.fits
    outhdul.flush()
    outhdul.close()
    hdu.close()

''''
import glob
from show_allred import show_allred,show_corr_total_flux

# basedir = '/allegro1/matisse/varga/'
basedir = '/data2/archive/data/MATISSE/'
reduced_dir = basedir + '/matisse_red5/'
targetdir = basedir +'/targets/'
targetlist_path = targetdir+'mat_target_list.txt'
target = 'HD_97048'
tpl_start = '2019-05-15T00_10_56'
'''
# target = 'HD_166191'
# tpl_start = '2019-09-19T00_39_18'
# target = 'HD_163296'
# tpl_start = '2019-05-06T08_19_51' #'2019-03-23T08_41_19' #
# target = 'RU_Lup'
# tpl_start = '2019-05-14T03_06_51'

#producing figures for the proposal (HD97048):
# band = 'L'
# target = 'HD_97048'
# tpl_start = '2019-05-15T00_10_56'
# input_fitsfiles = [targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_'+'coherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_'+'incoherent'+'_TARGET_FLUXCAL_INT_avg.fits']
# input_fitsfile_types = ['c','f']
# show_corr_total_flux(input_fitsfiles,input_fitsfile_types, outputdir=targetdir + '/' + target + '/', 
#     fn_pattern = '', verbose=False, save_png=True, save_eps=False,sel_wl=3.585,
#     bandwidth=0.03,file_type='',pro_catg='')

# #producing figures for the proposal (HD163296):
# band = 'L'
# target = 'HD_163296'
# tpl_start1 = '2019-03-23T08_41_19'
# tpl_start2 = '2019-05-06T08_19_51'
# input_fitsfiles = [targetdir + '/' + target + '/'+target+'_'+tpl_start1+'_'+band+'_'+'coherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# targetdir + '/' + target + '/'+target+'_'+tpl_start1+'_'+band+'_'+'incoherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# targetdir + '/' + target + '/'+target+'_'+tpl_start2+'_'+band+'_'+'coherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# #targetdir + '/' + target + '/'+target+'_'+tpl_start2+'_'+band+'_'+'incoherent'+'_TARGET_FLUXCAL_INT_avg.fits'
# ]
# input_fitsfile_types = ['c','f','c'] #,'f']
# show_corr_total_flux(input_fitsfiles,input_fitsfile_types, outputdir=targetdir + '/' + target + '/', 
#     fn_pattern = '', verbose=False, save_png=True, save_eps=False,sel_wl=3.6,
#     bandwidth=0.4,file_type='',pro_catg='')

# band = 'N'
# target = 'HD_163296'
# tpl_start1 = '2019-03-23T08_41_19'
# tpl_start2 = '2019-05-06T08_19_51'
# input_fitsfiles = [targetdir + '/' + target + '/'+target+'_'+tpl_start1+'_'+band+'_'+'coherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# targetdir + '/' + target + '/'+target+'_'+tpl_start1+'_'+band+'_'+'incoherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# targetdir + '/' + target + '/'+target+'_'+tpl_start2+'_'+band+'_'+'coherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# #targetdir + '/' + target + '/'+target+'_'+tpl_start2+'_'+band+'_'+'incoherent'+'_TARGET_FLUXCAL_INT_avg.fits'
# ]
# input_fitsfile_types = ['c','f','c'] #,'f']
# show_corr_total_flux(input_fitsfiles,input_fitsfile_types, outputdir=targetdir + '/' + target + '/', 
#     fn_pattern = '', verbose=False, save_png=True, save_eps=False,sel_wl=10.7,
#     bandwidth=0.4,file_type='',pro_catg='')

# band = 'L'
# target = 'HD_166191'
# tpl_start = '2019-09-19T00_39_18'
# input_fitsfiles = [targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_'+'coherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_'+'incoherent'+'_TARGET_FLUXCAL_INT_avg.fits']
# input_fitsfile_types = ['c','f']
# show_corr_total_flux(input_fitsfiles,input_fitsfile_types, outputdir=targetdir + '/' + target + '/', 
#     fn_pattern = '', verbose=False, save_png=True, save_eps=False,sel_wl=3.6,
#     bandwidth=0.4,file_type='',pro_catg='',annotate=True)

# band = 'N'
# target = 'HD_166191'
# tpl_start = '2019-09-19T00_39_18'
# input_fitsfiles = [targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_'+'coherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_'+'incoherent'+'_TARGET_FLUXCAL_INT_avg.fits']
# input_fitsfile_types = ['c','f']
# show_corr_total_flux(input_fitsfiles,input_fitsfile_types, outputdir=targetdir + '/' + target + '/', 
#     fn_pattern = '', verbose=False, save_png=True, save_eps=False,sel_wl=10.7,
#     bandwidth=0.4,file_type='',pro_catg='',annotate=True)

# #producing figures for the proposal (RU Lup):
# band = 'L'
# target = 'RU_Lup'
# tpl_start = '2019-05-14T03_06_51'
# input_fitsfiles = [targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_'+'coherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_'+'incoherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# ]
# input_fitsfile_types = ['c','f'] 
# show_corr_total_flux(input_fitsfiles,input_fitsfile_types, outputdir=targetdir + '/' + target + '/', 
#     fn_pattern = '', verbose=False, save_png=True, save_eps=False,sel_wl=3.6,
#     bandwidth=0.3,file_type='',pro_catg='')

# band = 'N'
# target = 'RU_Lup'
# tpl_start = '2019-05-14T03_06_51'
# input_fitsfiles = [targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_'+'coherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_'+'incoherent'+'_TARGET_FLUXCAL_INT_avg.fits',
# ]
# input_fitsfile_types = ['c','f'] 
# show_corr_total_flux(input_fitsfiles,input_fitsfile_types, outputdir=targetdir + '/' + target + '/', 
#     fn_pattern = '', verbose=False, save_png=True, save_eps=False,sel_wl=10.7,
#     bandwidth=0.4,file_type='',pro_catg='')


# bands = ['L'] #'N
# cohs = ['coherent'] #['coherent','incoherent']
# tag = 'TARGET_CAL' #'TARGET_FLUXCAL'
# do_calc_vis_from_corrflux = False
# # input_corrflux_file = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+'L'+'_coherent'+'_TARGET_COHVISCAL_INT_avg.fits'
# # outfile_path = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+'L'+'_coherent'+'_TARGET_COHVISCAL_INT_avg_3um6.fits'
# # # input_corrflux_file = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+'N'+'_coherent'+'_TARGET_COHVISCAL_INT_avg.fits'
# # # outfile_path = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+'N'+'_coherent'+'_TARGET_COHVISCAL_INT_avg_10um7.fits'
# # input_corrflux_file = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+'N'+'_incoherent'+'_TARGET_CAL_INT_avg.fits'
# # outfile_path = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+'N'+'_incoherent'+'_TARGET_CAL_INT_avg_8um5.fits'
# # input_corrflux_file = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+'L'+'_coherent'+'_TARGET_COHVISCAL_INT_avg.fits'
# # outfile_path = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+'L'+'_coherent'+'_TARGET_COHVISCAL_INT_avg_4um0.fits'
# # oifits_restrict_wavelengths(input_corrflux_file,outfile_path,4.,0.05) # 3.6,0.05) #10.7,0.1) #8.5,0.1
# # EXTERMINATE
# for band in bands:
#     print(band)
#     for coh in cohs:
#         infile_list = glob.glob(targetdir + '/' + target + '/*'+tpl_start+'_'+band+'_'+coh+'_'+tag+'_INT_0*.fits')
#         infile_list = sorted(infile_list)
#         outfile_path = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_'+coh+'_'+tag+'_INT_avg.fits'
#         # print(infile_list)
#         avg_oifits(infile_list,outfile_path)
#     if do_calc_vis_from_corrflux == True:
#         input_corrflux_file = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_coherent'+'_TARGET_FLUXCAL_INT_avg.fits'
#         input_totalflux_file = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_incoherent'+'_TARGET_FLUXCAL_INT_avg.fits'
#         outfile_path = targetdir + '/' + target + '/'+target+'_'+tpl_start+'_'+band+'_coherent'+'_TARGET_COHVISCAL_INT_avg.fits'
#         calc_vis_from_corrflux(input_corrflux_file,input_totalflux_file,outfile_path)
# inputdir = targetdir + '/' + target + '/'
# # print(inputdir)
# #show_allred(inputdir, outputdir=inputdir,fn_pattern='*'+tag+'_INT_avg', verbose=False,sel_wls=[3.2,3.535,3.585,3.6,4.0,10.7],bandwidths=[0.2,0.025,0.03,0.2,0.2,0.3])                      
# show_allred(inputdir, outputdir=inputdir+'/plots/',fn_pattern='*'+tag+'_INT_avg', verbose=False,fit_model=True,
# sel_wls=[3.53,3.58],bandwidths=[0.02,0.04])
# # sel_wls=[3.6,4.0,10.7],bandwidths=[0.2,0.2,0.3])
# if do_calc_vis_from_corrflux == True:
#     show_allred(inputdir, outputdir=inputdir+'/plots/',fn_pattern='*TARGET_COHVISCAL_INT_avg', verbose=False,sel_wls=[3.6,4.0,10.7],bandwidths=[0.2,0.2,0.3],fit_model=True)

# print('READY')
