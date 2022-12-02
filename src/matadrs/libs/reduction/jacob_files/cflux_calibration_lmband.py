__author__ = "Jacob Isbell"

import numpy as np
from astropy.io import fits
from sys import argv
import sys
import matplotlib.pyplot as plt
from scipy.special import jv, jve

script, sof, outdir = argv

est_from_cycle = True

def fopen(fname):
    hdu = fits.open(fname)
    return hdu

def read_spectrum(wl):
    f = open('../../calibrator_templates/HD120404_spectrum.txt','r')
    lines = f.readlines()
    wave, flux = [], []
    for l in lines:
        if l[0] != '#':
            data = l.split()
            wave.append(float(data[0]))
            flux.append(float(data[1]))

    myflux = np.interp(wl*1e6, wave, flux )
    return myflux


#first load the sof
targets = []
calibs = []
f = open(sof,'r')
lines = f.readlines()
for l in lines:
    if l[0] != '#':
        fname, ftype = l.split()
        if 'CALIB' in ftype:
            calibs.append( fopen(fname) )
        elif 'TARG' in ftype:
            targets.append( fopen(fname) )
        else:
            print('Unrecognized file in sof: %s'%(ftype))

if len(targets) != len(calibs):
    print('Number of calibs does not match number of targets. Exiting ...')
    sys.exit(1)


for i in range(len(targets)):
    t = targets[i]
    c = calibs[i]


    tpl_start = t[0].header['ESO TPL START']
    target = t[0].header['ESO OBS TARG NAME'].split('_')[1]
    detector = t[0].header['eso det name']

    bcd_config_t = t[0].header['HIERARCH ESO CFG BCD MODE'].lower()
    bcd_config_c = c[0].header['HIERARCH ESO CFG BCD MODE'].lower()

    if bcd_config_c != bcd_config_t:
        print('Please order BCDs the same for targ and calib. Exiting ...')
        sys.exit(1)

    keys = [['OI_VIS2', 'vis2data', 'vis2err'], ['OI_VIS', 'visamp', 'visamperr'] ]

    for k in keys[:]:
        #load the data in the proper fits extension
        ext, datakey, errkey = k
        print(ext, datakey, errkey)

        targ_val = t[ext].data[datakey]
        targ_err = t[ext].data[errkey]

        #est_from_cycle calculates the statistical errors on vis from the variations in the cycles of the obs
        if est_from_cycle and len(targ_val)>6:
            groups = [[],[],[],[],[],[]]
            for b in range(len(targ_val)):
                groups[b%6].append(targ_val[b])

            for b in range(len(targ_err)):
                targ_err[b] = np.std(groups[b%6],0)

        #calculate relative error for later
        targ_relerr = targ_err / targ_val

        #load the calibrator data
        calib_val = c[ext].data[datakey]
        print(t[0].header['eso obs targ name'],c[0].header['eso obs targ name'], c[0].header['ra'],c[0].header['dec'], t[0].header['ra'],t[0].header['dec'])
        calib_err = c[ext].data[datakey]
        calib_relerr = calib_err / calib_val

        u = c[ext].data['ucoord']
        v = c[ext].data['vcoord']

        #if there are more cycles of targ than calib, simply repeat the calib ones
        if len(targ_val) > len(calib_val):
            nrep = len(targ_val) // len(calib_val)
            new_val = np.copy(calib_val)
            new_err = np.copy(calib_err)
            new_u = np.copy(u)
            new_v = np.copy(v)
            for j in range(nrep-1):
                new_val = np.append(new_val, calib_val)
                new_err = np.append(new_err, calib_err)
                new_u = np.append(new_u, u)
                new_v = np.append(new_v, v)
            new_relerr = new_err / new_val

            calib_val = new_val.reshape(targ_val.shape)
            calib_relerr = new_relerr.reshape(targ_val.shape)
            calib_err = new_err.reshape(targ_val.shape)
            u = new_u.reshape(t[ext].data['ucoord'].shape )
            v = new_v.reshape(t[ext].data['vcoord'].shape )

        #need a flux template.... for now assume hot blackbody with no spectral features
        wl = t[3].data['eff_wave']
        rayleigh_jeans_approx = np.power(wl, -2)
        ######################
        ## set this value! ###
        ######################
        flux_scaling = 2.28 #Jy -- L band flux of your source
        cal_spectrum = flux_scaling * rayleigh_jeans_approx / np.mean(rayleigh_jeans_approx)


        #if you have a known stellar spectrum, use it here
        #cal_spectrum = read_spectrum(wl)   #in Jy,  from Roy van Boekel

        if( len(calib_val[0]) != len(targ_val[0]) ):
            print('wtf')
            print(wl)
            print(c['oi_wavelength'].data['eff_wave']  )
            print(np.interp(wl, c['oi_wavelength'].data['eff_wave'], calib_val[0]))
        temp_calib_val = []
        temp_calib_err = []#calib_err.copy()
        for b in range(len(calib_val)):
            temp_calib_val.append( np.interp(wl[::-1], c['oi_wavelength'].data['eff_wave'][::-1], calib_val[b][::-1] )[::-1] )
            temp_calib_err.append( np.interp(wl[::-1], c['oi_wavelength'].data['eff_wave'][::-1], calib_err[b][::-1] )[::-1] )

        calib_val=np.array(temp_calib_val)
        calib_err=np.array(temp_calib_err)
        calib_relerr = calib_err / calib_val

        #need to correct for resolution effects!!!!!
        ######################
        ## set this value! ###
        ######################
        d = 2.353381 #mas, diameter of calibrator star
        bls = np.sqrt( np.square(u) + np.square(v)  )
        #use uniform disk Fourier Transform to correct the visibilities
        for b in range(len(calib_val)):
            vis = np.absolute( 2 * jv(1, np.pi*d/1000/206265* bls[b]/wl  ) / (np.pi*d/1000/206265* bls[b]/wl)   )
            calib_val[b] /= vis

        cal_data = targ_val / calib_val * cal_spectrum#Jy
        cal_relerr = np.sqrt( np.square(targ_relerr))# + np.square(calib_relerr) ) #if uncertainty of spectrum is known, you can uncomment this
        cal_err = cal_data * calib_relerr

        #quick view plots
        if ext == "OI_VIS":
            fig = plt.figure()
            for b in range(len(targ_val)):
                #plt.errorbar(wl, targ_val[b],yerr=targ_err[b], c='blue')
                #plt.plot(wl, calib_val[b], c='red')
                plt.errorbar(wl*1e6, cal_data[b],yerr=cal_err[b], c='k')
                #plt.plot(wl, cal_spectrum,'r--' )
                s = np.where(np.logical_and( wl >3.5e-6, wl<3.7e-6  ))[0]
            plt.ylim([0,5*np.mean(cal_data[0][s])      ])
            plt.xlabel('wl [micron]')
            plt.ylabel('flux [Jy]')
            s = np.where(np.logical_and( wl >3.45e-6, wl<3.55e-6  ))[0]
            plt.title(r'Flux at 3.5um: $%.3f \pm %.3f$ Jy'%( np.mean(cal_data[:,s] ), np.std(cal_data[:,s] )  )  )
            #plt.show()
            outname = outdir + '/' + target.lower()+'_'+tpl_start+'_' +detector.upper() + '_' + bcd_config_t + '.png'
            plt.savefig(outname, bbinches='tight')
            plt.close()

        t[ext].data[datakey] = cal_data
        t[ext].data[errkey] = cal_err

    outname = outdir + '/' + target.lower()+'_'+tpl_start+'_' +detector.upper() + '_' + bcd_config_t + '.fits'
    print(outname)
    t.writeto(outname, overwrite=True)
