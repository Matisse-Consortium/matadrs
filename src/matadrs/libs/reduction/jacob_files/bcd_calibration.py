#!/usr/bin/env python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from astropy.io import fits
from glob import glob
from pathlib import Path
from sys import argv

mpl.rcParams['font.family']='serif'
mpl.rcParams['xtick.direction']='in'
mpl.rcParams['xtick.top']=True
mpl.rcParams['ytick.direction']='in'
mpl.rcParams['ytick.right']=True

def compare_cphase_OO(fname: Path) -> None:
    hdu = fits.open(fname)

    cphase = hdu['oi_t3'].data['t3phi']
    cphase_err = hdu['oi_t3'].data['t3phierr']
    wl = hdu['oi_wavelength'].data['eff_wave']

    #print(hdu['oi_t3'].data['sta_index'] - 31)
    """
    [[2 3 4]
     [1 2 3]
     [1 2 4]
     [1 3 4]]

    """

    test_cphase = cphase[3] - cphase[2] + cphase[1]
    test_cphase_err = cphase_err[1] + cphase_err[3] - cphase_err[2]

    fig, axarr = plt.subplots(2,2)
    ax,bx,cx,dx = axarr.flatten()

    ax.errorbar(wl, cphase[0], yerr=cphase_err[0] )
    ax.errorbar(wl, test_cphase, yerr=test_cphase_err)
    #ax.errorbar(wl, np.array(-test_cphase+cphase[0]) )

    test_cphase = cphase[0] + cphase[3] + cphase[2]
    test_cphase_err = cphase_err[0] - cphase_err[3] + cphase_err[2]
    bx.errorbar(wl, cphase[1], yerr=cphase_err[1] )
    bx.errorbar(wl, test_cphase, yerr=test_cphase_err)
    #bx.errorbar(wl, np.array(test_cphase-cphase[1]) )

    test_cphase = cphase[1] + cphase[3] + cphase[0]
    test_cphase_err = cphase_err[1] + cphase_err[3] - cphase_err[0]
    cx.errorbar(wl, cphase[2], yerr=cphase_err[2] )
    cx.errorbar(wl, test_cphase, yerr=test_cphase_err)

    test_cphase = cphase[0] + cphase[2] + cphase[1]
    test_cphase_err = cphase_err[0] + cphase_err[2] - cphase_err[1]
    dx.errorbar(wl, cphase[3], yerr=cphase_err[3] )
    dx.errorbar(wl, test_cphase, yerr=test_cphase_err)

    plt.show()

def spectral_binning(arr: np.array, binwidth: int = 1, phase: bool = False, err: bool = False) -> np.array:
    new_arr = []
    new_err = []
    #print(arr.shape)
    if arr.ndim > 1:
        for subarr in arr:
            temp = []
            temperr = []
            for i in range(0,len(subarr)-binwidth,binwidth):
                if phase:
                    phase_arr = np.exp(1j*np.radians(subarr[i:i+binwidth]) )
                    if err:
                        temp.append(np.std(np.angle(phase_arr,deg=True) ))
                    else:
                        temp.append(np.angle(np.nanmedian(phase_arr),deg=True ))
                else:
                    temp.append(np.nanmedian(subarr[i:i+binwidth] ) )
                    temperr.append(np.nanstd(subarr[i:i+binwidth] ))
                    print(temperr)
            new_arr.append(temp)
            new_err.append(temperr)
    else:
        for i in range(0,len(arr)-binwidth,binwidth):
            new_arr.append(np.nanmedian(arr[i:i+binwidth] ) )
            new_err.append(np.std(arr[i:i+binwidth] ) )
            #print(np.std(arr[i:i+binwidth] ))
    return np.array(new_arr)#, np.array(new_err)


def cphase_calib(f1: Path, f2: Path, f3: Path, f4: Path) -> None:
    hdu1 = fits.open(f1)

    cphase1 = hdu1['oi_t3'].data['t3phi']
    # print(cphase1.shape, 'testing')
    cphase1_err = hdu1['oi_t3'].data['t3phierr']
    wl = hdu1['oi_wavelength'].data['eff_wave']
    loops = hdu1['OI_T3'].data['sta_index']
    tel_names = hdu1[2].data['tel_name']
    sta_name = hdu1[2].data['sta_index']
    all_tels = ['A0', 'B2', 'C0', 'D1'] + ['K0', 'G1', 'D0', 'J3'] + ['A0', 'G1', 'J2', 'J3'] + ['UT1', 'UT2', 'UT3', 'UT4']    # Different baseline-configurations short-, medium-, large AT, UT
    all_stas = [1,  5, 13, 10] + [28, 18, 13, 24] + [1, 18, 23, 24] + [32, 33, 34, 35]                                          # 'sta_index'of telescopes
    telescopes = []

    hdu2 = fits.open(f2)
    cphase2 = hdu2['oi_t3'].data['t3phi']
    cphase2_err = hdu2['oi_t3'].data['t3phierr']

    hdu3 = fits.open(f3)
    cphase3 = hdu3['oi_t3'].data['t3phi']
    cphase3_err = hdu3['oi_t3'].data['t3phierr']

    hdu4 = fits.open(f4)
    cphase4 = hdu4['oi_t3'].data['t3phi']
    cphase4_err = hdu4['oi_t3'].data['t3phierr']



    """
    outout
    [[2 3 4]
     [1 2 3]
     [1 2 4]
     [1 3 4]]
    inin
    [[1 4 3]
     [2 1 4]
     [2 1 3]
     [2 4 3]]
    inout
    [[2 4 3]
     [1 2 4]
     [1 2 3]
     [1 4 3]]
    outin
    [[1 3 4]
     [2 1 3]
     [2 1 4]
    [2 3 4]]

    """

    corr_cphases = []
    corr_cphaseerr = []

    #print(cphase1.shape[0] // 4, cphase2.shape[0] // 4)
    #length = np.min([cphase1.shape[0],cphase2.shape[0]]) //4
    s = np.argmin([cphase1.shape[0],cphase2.shape[0]])
    length = [cphase1.shape[0],cphase2.shape[0]][s] //4
    for kk in range(length):
        index1 = 0 + kk*4
        index2 = 3 + kk*4
        index3 = 0 + kk*4
        index4 = 3 + kk*4

        #ax.errorbar(wl*1e6, -cphase3[index3],yerr=cphase3_err[index3], alpha=0.25, color='grey', zorder=-1)
        #ax.errorbar(wl*1e6, cphase4[index4],yerr=cphase4_err[index4], alpha=0.25, color='grey', zorder=-1)
        mean_val =  np.nanmean( [cphase1[index1],-cphase2[index2],-cphase3[index3],cphase4[index4] ],0)
        std_val =  np.nanstd( [cphase1[index1],-cphase2[index2], -cphase3[index3],cphase4[index4]  ],0)
        #print(np.nanmean(cphase1[index1] ), np.nanmean(cphase2[index2] ) )
        #mean_val =  np.nanmean( [cphase1[index1],-cphase2[index2] ],0)
        #mean_val =  (cphase1[index1] + cphase2[index2]) / 2
        #std_val =  np.nanstd( [cphase1[index1],-cphase2[index2] ],0)

        corr_cphases.append(mean_val)
        corr_cphaseerr.append(std_val)

        index1 = 1 + kk*4
        index2 = 2 + kk*4
        index3 = 2 + kk*4
        index4 = 1 + kk*4
        mean_val =  np.nanmean( [cphase1[index1],-cphase2[index2],cphase3[index3],-cphase4[index4] ],0)
        std_val =  np.nanstd( [cphase1[index1],-cphase2[index2], cphase3[index3],-cphase4[index4]  ],0)
        #mean_val =  np.nanmean( [cphase1[index1],-cphase2[index2] ],0)
        #std_val =  np.nanstd( [cphase1[index1],-cphase2[index2] ],0)
        corr_cphases.append(mean_val)
        corr_cphaseerr.append(std_val)

        index1 = 2 + kk*4
        index2 = 1 + kk*4
        index3 = 1 + kk*4
        index4 = 2 + kk*4
        mean_val =  np.nanmean( [cphase1[index1],-cphase2[index2],cphase3[index3],-cphase4[index4] ],0)
        std_val =  np.nanstd( [cphase1[index1],-cphase2[index2], cphase3[index3],-cphase4[index4]  ],0)
        #mean_val =  np.nanmean( [cphase1[index1],-cphase2[index2] ],0)
        #std_val =  np.nanstd( [cphase1[index1],-cphase2[index2] ],0)
        corr_cphases.append(mean_val)
        corr_cphaseerr.append(std_val)

        index1 = 3 + kk*4
        index2 = 0 + kk*4
        index3 = 0 + kk*4
        index4 = 3 + kk*4
        mean_val =  np.nanmean( [cphase1[index1],-cphase2[index2],-cphase3[index3],cphase4[index4] ],0)
        std_val =  np.nanstd( [cphase1[index1],-cphase2[index2], -cphase3[index3],cphase4[index4]  ],0)
        #mean_val =  np.nanmean( [cphase1[index1],-cphase2[index2] ],0)
        #std_val =  np.nanstd( [cphase1[index1],-cphase2[index2] ],0)

        corr_cphases.append(mean_val)
        corr_cphaseerr.append(std_val)


    corr_cphases = np.array(corr_cphases)
    corr_cphaseerr = np.array(corr_cphaseerr)

    #print(hdu1[0].header['eso tpl start'], corr_cphases, corr_cphaseerr )
    #print(corr_cphaseerr.shape)
    #hdus = [hdu1,hdu2]
    hdus = [hdu1,hdu2,hdu3,hdu4]

    files = [f1,f2]

    try:
        myhdu = fits.open(files[s].split('.fits')[0].replace('out-out','').replace('in-in','').replace('out-in','').replace('in-out','') + 'bcd_calibratedTEST.fits')
    except:
        myhdu = hdus[s]



    myhdu['oi_t3'].data['t3phierr'] = corr_cphaseerr
    myhdu['oi_t3'].data['t3phi'] = corr_cphases
    myhdu.writeto(files[s].split('.fits')[0].replace('out-out','').replace('in-in','').replace('out-in','').replace('in-out','') + 'bcd_calibratedTEST.fits', overwrite=True )

    '''
    fig,axarr = plt.subplots(2,2)
    groups = [[],[],[],[]]
    for i in range(len(corr_cphases)):
        ax = axarr.flatten()[i%4]
        ax.errorbar(wl*1e6, corr_cphases[i], yerr=corr_cphaseerr[i],zorder=1,alpha=0.5)
        groups[i%4].append(corr_cphases[i])
        ax.set_ylim([-181,181])

    for i,g in enumerate(groups):
        ax = axarr.flatten()[i]
        print(np.nanmean(g,0) )
        ax.errorbar(wl*1e6, np.nanmean(g,0), yerr=np.nanstd(g,0),ls='--', lw=2,zorder=2,color='k' )

    #plt.show()
    #plt.close()

    # Creates folder for plots if it does not exist
    if not os.path.isdir("../plots/bcd/"):
        os.makedirs("../plots/bcd/")

    plt.savefig('../plots/bcd/%s.png'%(hdu1[0].header['eso tpl start'].replace(":","_") ))
    # plt.show()
    plt.close()
    '''
    return


def vis2_calib(f1: Path, f2: Path, f3: Path, f4: Path, cflux: bool = False) -> None:
    hdu1 = fits.open(f1)

    key = 'oi_vis2'
    key_data = 'vis2data'
    key_err = 'vis2err'
    if cflux:
        key = 'oi_vis'
        key_data = 'visamp'
        key_err = 'visamperr'

    vis21 = hdu1[key].data[key_data]
    vis21_err = hdu1[key].data[key_err]
    wl = hdu1['oi_wavelength'].data['eff_wave']
    loops = hdu1['OI_vis2'].data['sta_index']
    tel_names = hdu1[2].data['tel_name']
    sta_name = hdu1[2].data['sta_index']
    all_tels = ['A0', 'B2', 'C0', 'D1'] + ['K0', 'G1', 'D0', 'J3'] + ['A0', 'G1', 'J2', 'J3'] + ['UT1', 'UT2', 'UT3', 'UT4']    # Different baseline-configurations short-, medium-, large AT, UT
    all_stas = [1,  5, 13, 10] + [28, 18, 13, 24] + [1, 18, 23, 24] + [32, 33, 34, 35]                                          # 'sta_index'of telescopes
    telescopes = []

    hdu2 = fits.open(f2)
    #vis22 = hdu2['oi_vis2'].data['vis2data']
    #vis22_err = hdu2['oi_vis2'].data['vis2err']
    vis22 = hdu2[key].data[key_data]
    vis22_err = hdu2[key].data[key_err]

    hdu3 = fits.open(f3)
    #vis23 = hdu3['oi_vis2'].data['vis2data']
    #vis23_err = hdu3['oi_vis2'].data['vis2err']
    vis23 = hdu3[key].data[key_data]
    vis23_err = hdu3[key].data[key_err]

    hdu4 = fits.open(f4)
    #vis24 = hdu4['oi_vis2'].data['vis2data']
    #vis24_err = hdu4['oi_vis2'].data['vis2err']
    vis24 = hdu4[key].data[key_data]
    vis24_err = hdu4[key].data[key_err]



    """
    outout
    [[34, 35],
   [32, 33],
   [33, 34],
   [33, 35],
   [32, 34],
   [32, 35]]
    inin
    [35, 34],
   [33, 32],
   [32, 35],
   [32, 34],
   [33, 35],
   [33, 34]]
   outin
   [34, 35],
   [33, 32],
   [32, 34],
   [32, 35],
   [33, 34],
   [33, 35]
   inout
   [35, 34],
   [32, 33],
   [33, 35],
   [33, 34],
   [32, 35],
   [32, 34]

    """

    corr_vis2s = []
    corr_vis2err = []
    print(vis21.shape[0] // 6, 'vis2')
    s = np.argmin([vis21.shape[0],vis22.shape[0]])
    length = [vis21.shape[0], vis22.shape[0]][s] //6
    for kk in range(length):
        index1 = 0 + kk*6
        index2 = 0 + kk*6
        index3 = 0 + kk*6
        index4 = 0 + kk*6

        #ax.errorbar(wl*1e6, -vis23[index3],yerr=vis23_err[index3], alpha=0.25, color='grey', zorder=-1)
        #ax.errorbar(wl*1e6, vis24[index4],yerr=vis24_err[index4], alpha=0.25, color='grey', zorder=-1)
        mean_val =  np.nanmedian( [vis21[index1],vis22[index2],vis23[index3],vis24[index4] ],0)
        std_val =  np.nanstd( [vis21[index1],vis22[index2], vis23[index3],vis24[index4]  ],0)
        #mean_val =  np.nanmedian( [vis21[index1],vis22[index2] ],0)
        #std_val =  np.nanstd( [vis21[index1],vis22[index2] ],0)
        corr_vis2s.append(mean_val)
        corr_vis2err.append(std_val)

        index1 = 1 + kk*6
        index2 = 1 + kk*6
        index3 = 2 + kk*6
        index4 = 1 + kk*6
        mean_val =  np.nanmedian( [vis21[index1],vis22[index2],vis23[index3],vis24[index4] ],0)
        std_val =  np.nanstd( [vis21[index1],vis22[index2], vis23[index3],vis24[index4]  ],0)
        #mean_val =  np.nanmedian( [vis21[index1],vis22[index2] ],0)
        #std_val =  np.nanstd( [vis21[index1],vis22[index2] ],0)
        corr_vis2s.append(mean_val)
        corr_vis2err.append(std_val)

        index1 = 2 + kk*6
        index2 = 5 + kk*6
        index3 = 4 + kk*6
        index4 = 3 + kk*6
        #mean_val =  np.nanmedian( [vis21[index1],-vis22[index2],vis23[index3],-vis24[index4] ],0)
        #std_val =  np.nanstd( [vis21[index1],-vis22[index2], vis23[index3],-vis24[index4]  ],0)
        #mean_val =  np.nanmedian( [vis21[index1],vis22[index2] ],0)
        #std_val =  np.nanstd( [vis21[index1],vis22[index2] ],0)
        mean_val =  np.nanmedian( [vis21[index1],vis22[index2],vis23[index3],vis24[index4] ],0)
        std_val =  np.nanstd( [vis21[index1],vis22[index2], vis23[index3],vis24[index4]  ],0)
        corr_vis2s.append(mean_val)
        corr_vis2err.append(std_val)

        index1 = 3 + kk*6
        index2 = 4 + kk*6
        index3 = 5 + kk*6
        index4 = 2 + kk*6
        #mean_val =  np.nanmedian( [vis21[index1],-vis22[index2],-vis23[index3],vis24[index4] ],0)
        #std_val =  np.nanstd( [vis21[index1],-vis22[index2], -vis23[index3],vis24[index4]  ],0)
        #mean_val =  np.nanmedian( [vis21[index1],vis22[index2] ],0)
        #std_val =  np.nanstd( [vis21[index1],vis22[index2] ],0)
        mean_val =  np.nanmedian( [vis21[index1],vis22[index2],vis23[index3],vis24[index4] ],0)
        std_val =  np.nanstd( [vis21[index1],vis22[index2], vis23[index3],vis24[index4]  ],0)
        corr_vis2s.append(mean_val)
        corr_vis2err.append(std_val)

        index1 = 4 + kk*6
        index2 = 3 + kk*6
        index3 = 2 + kk*6
        index4 = 5 + kk*6
        #mean_val =  np.nanmedian( [vis21[index1],-vis22[index2],-vis23[index3],vis24[index4] ],0)
        #std_val =  np.nanstd( [vis21[index1],-vis22[index2], -vis23[index3],vis24[index4]  ],0)
        mean_val =  np.nanmedian( [vis21[index1],vis22[index2],vis23[index3],vis24[index4] ],0)
        std_val =  np.nanstd( [vis21[index1],vis22[index2], vis23[index3],vis24[index4]  ],0)
        #mean_val =  np.nanmedian( [vis21[index1],vis22[index2] ],0)
        #std_val =  np.nanstd( [vis21[index1],vis22[index2] ],0)
        corr_vis2s.append(mean_val)
        corr_vis2err.append(std_val)

        index1 = 5 + kk*6
        index2 = 2 + kk*6
        index3 = 3 + kk*6
        index4 = 4 + kk*6
        #mean_val =  np.nanmedian( [vis21[index1],-vis22[index2],-vis23[index3],vis24[index4] ],0)
        #std_val =  np.nanstd( [vis21[index1],-vis22[index2], -vis23[index3],vis24[index4]  ],0)
        mean_val =  np.nanmedian( [vis21[index1],vis22[index2],vis23[index3],vis24[index4] ],0)
        std_val =  np.nanstd( [vis21[index1],vis22[index2], vis23[index3],vis24[index4]  ],0)
        #mean_val =  np.nanmedian( [vis21[index1],vis22[index2] ],0)
        #std_val =  np.nanstd( [vis21[index1],vis22[index2] ],0)
        corr_vis2s.append(mean_val)
        corr_vis2err.append(std_val)

    corr_vis2s = np.array(corr_vis2s)
    corr_vis2err = np.array(corr_vis2err)

    #print(hdu1[0].header['eso tpl start'], corr_vis2s, corr_vis2err )
    print(corr_vis2err.shape)

    hdus = [hdu1,hdu2]

    files = [f1,f2]

    try:
        myhdu = fits.open(files[s].split('.fits')[0].replace('out-out','').replace('in-in','').replace('out-in','').replace('in-out','') + 'bcd_calibratedTEST.fits')
    except:
        myhdu = hdus[s]



    myhdu[key].data[key_data] = corr_vis2s
    myhdu[key].data[key_err] = corr_vis2err
    #if cflux:
    #	myhdu['oi_vis2'].data['vis2data'] = corr_vis2s
    #	myhdu['oi_vis2'].data['vis2err'] = corr_vis2err


    myhdu.writeto(files[s].split('.fits')[0].replace('out-out','').replace('in-in','').replace('out-in','').replace('in-out','') + 'bcd_calibratedTEST.fits', overwrite=True )

    '''
    fig,axarr = plt.subplots(2, 3)

    groups = [[],[],[],[],[],[]]
    for i, o in enumerate(corr_vis2s):
        ax = axarr.flatten()[i%6]
        ax.errorbar(wl*1e6, o, yerr=corr_vis2err[i], zorder=1, alpha=0.5)
        ax.set_ylim([0, 0.2])
        groups[i%6].append(o)
        if cflux:
            ax.set_ylim([0, 2])

    for j, g in enumerate(groups):
        ax = axarr.flatten()[j]
        print(np.nanmedian(g, 0) )
        ax.errorbar(wl*1e6, np.nanmedian(g, 0), yerr=np.nanstd(g, 0),ls='--', lw=2, zorder=2, color='k')

    # plt.show()
    plt.close()

    # Creates folder for plots if it does not exist
    if not os.path.isdir("../plots/bcd/"):
        os.makedirs("../plots/bcd/")

    plt.savefig('../plots/bcd/%s.png'%(hdu1[0].header['eso tpl start'].replace(":","_") ))
    # plt.show()
    plt.close()
    '''
    return


def main(fdir: Path) -> None:
    """
    files_oo = np.sort(glob(fdir+'/*out-out.fits'))
    files_ii = np.sort(glob(fdir+'/*in-in.fits'))
    files_io = np.sort(glob(fdir+'/*in-out.fits'))
    files_oi = np.sort(glob(fdir+'/*out-in.fits'))
    """
    files_oo = np.sort(glob(fdir+'/*0001.fits'))
    files_ii = np.sort(glob(fdir+'/*0002.fits'))
    files_io = np.sort(glob(fdir+'/*0003.fits'))
    files_oi = np.sort(glob(fdir+'/*0004.fits'))

    for i, o in enumerate(files_oo[:]):
        cphase_calib(o, files_ii[i], files_io[i], files_oi[i])
        vis2_calib(o, files_ii[i], files_oo[i], files_oo[i], cflux=False)


if __name__ == "__main__":
    try:
        script, fdir = argv
    except:
        print("Please enter the directory containing the calibrated data. E.g. 'python3 bcd_calibration.py /path/to/directory/'")
        sys.exit(1)

    main(fdir)

