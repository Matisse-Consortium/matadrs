from mat_tools import mat_autoPipeline as mat_autoPipeline
import matplotlib
matplotlib.use('Agg')
import os
import shutil
import glob
import numpy as np
#from mcdb import Mcdb as mc
# from mcdb import VioletasExp as ve
# from mcdb import matutil as mu
# from mcdb import wutil as wu
from show_allred import show_allred, show_allred_mosaic
from pk2oifits import pk_2_oifits
print(mat_autoPipeline.__file__)

DATADIR = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/GTO/hd142666/RAW"
RESDIR_MAIN = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/GTO/hd142666/PRODUCTS/20190514"
PROGRAMDIR = '/allegro6/matisse/varga/pro/'

#test_red/matisse_red7: N spectralbinning_N = '7'
#test_red/matisse_red7_test: N spectralbinning_N = '147'

do_L = False
do_N = True
do_reduction = True
do_plot = True

pipeline = 'drs' #'drs', 'ews' #

spectralbinning_L = '5' #DRS option (tried value = 5 and 10 - matisse_redM), default for LR: '1'
spectralbinning_N = '7' #21 #DRS option, default for LR: '7' (tried value = 49 - matisse_redM and some matisse_red6)

ews_modes = ['ews']

# incoherent processing (default)
# coherent processing: corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2/
drs_modes = [
    'coherent',
    'incoherent'
    ]

if pipeline == 'drs':
    modes = drs_modes
if pipeline == 'ews':
    modes = ews_modes

paramL_list = [
'/spectralBinning='+spectralbinning_L+'/corrFlux=TRUE/useOpdMod=FALSE/coherentAlgo=2/compensate="[pb,rb,nl,if,bp,od]"/', #for coherent visibilities: '/corrFlux=FALSE/useOpdMod=FALSE/coherentAlgo=2/compensate="[pb,rb,nl,if,bp,od]"/'
    '/spectralBinning='+spectralbinning_L+'/compensate="pb,rb,nl,if,bp,od"/'
    ]

paramN_list_AT = [
'/replaceTel=3/corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2/spectralBinning='+spectralbinning_N,
      '/replaceTel=3/spectralBinning='+spectralbinning_N
    ]
paramN_list_UT = [
'/replaceTel=0/corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2/spectralBinning='+spectralbinning_N,
      '/replaceTel=0/spectralBinning='+spectralbinning_N
    ]


obs_lst=[
{'night':'2019-03-22', 'tpl_start':'2019-03-23T07:44:51', 'tel':'ATs','diL':'LOW','diN':'LOW'},  #IRAS_13481-6124, SCI_IRAS13481-6124
{'night':'2019-03-22', 'tpl_start':'2019-03-23T08:15:36', 'tel':'ATs','diL':'LOW','diN':'LOW'},  #zet_Ara, CAL_IRAS13481-6124_zet_Ara
]

N_obs = len(obs_lst)

if do_L:
    skip_L = 0
else:
    skip_L = 1
if do_N:
    skip_N = 0
else:
    skip_N = 1

for i in range(N_obs):
    TPL_START=obs_lst[i]['tpl_start']
    NIGHT=obs_lst[i]['night']
    TEL=obs_lst[i]['tel']
    DIL=obs_lst[i]['diL']
    DIN=obs_lst[i]['diN']
    # ----------define the folders for saving data-----------------
    RAWDIR = DATADIR + '/' + NIGHT
    CALIBDIR = '/allegro6/matisse/varga/CalibMap'
    CALIBDIR2 = '/allegro6/matisse/drs/calibmaps/CalibMapMayJuneJuly2019'

    # ----------run the pipeline-------------------------------
    if pipeline == 'drs':
        for j in range(len(modes)):
            RESDIR = RESDIR_MAIN+'/'+modes[j]+'/' + NIGHT + '/' + TPL_START.replace(':', '_')
            if not os.path.exists(RESDIR):
                os.makedirs(RESDIR)
            if do_reduction:
                if do_L:
                    RESDIR_L = glob.glob(RESDIR + '/Iter1/*HAWAII*/')
                    if RESDIR_L:
                        shutil.rmtree(RESDIR_L[0])
                    soffiles = glob.glob(RESDIR + '/Iter1/*HAWAII*.sof*')
                    if soffiles:
                        for file in soffiles:
                            os.remove(file)
                if do_N:
                    RESDIR_N = glob.glob(RESDIR + '/Iter1/*AQUARIUS*/')
                    if RESDIR_N:
                        shutil.rmtree(RESDIR_N[0])
                    soffiles = glob.glob(RESDIR + '/Iter1/*AQUARIUS*.sof*')
                    if soffiles:
                        for file in soffiles:
                            os.remove(file)

                if TEL == 'UTs':
                    paramN_list = paramN_list_UT
                if TEL == 'ATs':
                    paramN_list = paramN_list_AT
                if 'HIGH' in DIN:
                    for ih in range(len(paramN_list)):
                        paramN_list[ih] = paramN_list[ih].replace('spectralBinning=7','spectralBinning=49') 
                else: 
                    for ih in range(len(paramN_list)):
                        paramN_list[ih] = paramN_list[ih].replace('spectralBinning=49','spectralBinning=7') 
                # first try to find calibration files in RAWDIR+'/calibration_files/'
                if os.path.exists(RAWDIR+'/calibration_files/'):
                    res = mat_autoPipeline.mat_autoPipeline(dirRaw=RAWDIR, dirResult=RESDIR, dirCalib=RAWDIR+'/calibration_files/', nbCore=6, tplstartsel=TPL_START,
                                                  resol='', paramL=paramL_list[j], paramN=paramN_list[j], overwrite=0, maxIter=1,
                                                  skipL = skip_L, skipN = skip_N)
                    if res == 2: #if missing calibration
                        # if calibration files were not found, then use general calibration directory (CALIBDIR)
                        res = mat_autoPipeline.mat_autoPipeline(dirRaw=RAWDIR, dirResult=RESDIR, dirCalib=CALIBDIR, nbCore=6, tplstartsel=TPL_START,
                                                  resol='', paramL=paramL_list[j], paramN=paramN_list[j], overwrite=0, maxIter=1,
                                                  skipL = skip_L, skipN = skip_N)
                else:
                    # if there is no calibration directory within the night folder, then use general calibration directory (CALIBDIR)
                    res = mat_autoPipeline.mat_autoPipeline(dirRaw=RAWDIR, dirResult=RESDIR, dirCalib=CALIBDIR, nbCore=6, tplstartsel=TPL_START,
                                                  resol='', paramL=paramL_list[j], paramN=paramN_list[j], overwrite=0, maxIter=1,
                                                  skipL = skip_L, skipN = skip_N)
                    if res == 2:
                        # if calibration files were not found, then use Walter's calibration directory (CALIBDIR2)
                        res = mat_autoPipeline.mat_autoPipeline(dirRaw=RAWDIR, dirResult=RESDIR, dirCalib=CALIBDIR2, nbCore=6, tplstartsel=TPL_START,
                                                  resol='', paramL=paramL_list[j], paramN=paramN_list[j], overwrite=0, maxIter=1,
                                                  skipL = skip_L, skipN = skip_N)

    if pipeline == 'ews':
        for j in range(len(modes)):
            RESDIR = RESDIR_MAIN+'/'+modes[j]+'/' + NIGHT + '/' + TPL_START.replace(':', '_')
            if not os.path.exists(RESDIR):
                os.makedirs(RESDIR)
            current_dir = os.getcwd()
            os.chdir(RAWDIR)
            if not os.path.exists(RAWDIR+'/Tpl'):
                try:
                    mu.createTplDirectories(mother='Tpl')
                except FileExistsError as e:
                    print(e)

            if do_reduction:
                if do_L:
                    # create a directory for the reduced data
                    RESDIR_BAND = RESDIR + '/Iter1/mat_raw_estimates.'+TPL_START.replace(':', '_')+'.HAWAII-2RG.rb/'
                    if os.path.exists(RESDIR_BAND):
                        shutil.rmtree(RESDIR_BAND)
                    os.makedirs(RESDIR_BAND)
                    os.chdir(RAWDIR+'/Tpl')

                    try: 
                        try:
                            rdata = ve.processTemplate(TPL_START+'.tpl',band='LM',dophot=False)
                            if(rdata is not None):
                                ofile = RESDIR_BAND+'/'+TPL_START.replace(':', '_')+'.tpl.pk'
                                wu.msave(ofile, rdata)
                                pk_2_oifits(RESDIR_BAND+'/'+TPL_START.replace(':', '_')+'.tpl.pk', RESDIR_BAND,oifits_template=PROGRAMDIR+'/template_oifits.fits')
                        except FileNotFoundError as e:
                            print(e)
                    except IndexError as e:
                        print(e)

                if do_N:
                    # create a directory for the reduced data
                    RESDIR_BAND = RESDIR + '/Iter1/mat_raw_estimates.'+TPL_START.replace(':', '_')+'.AQUARIUS.rb/'
                    
                    if os.path.exists(RESDIR_BAND):
                        shutil.rmtree(RESDIR_BAND)
                    os.makedirs(RESDIR_BAND)
                    
                    os.chdir(RAWDIR+'/Tpl')

                    if 'HIGH' in DIN:
                        #ws = 21
                        #ws = 11 #!new value (from Jul 2021)!
                        #ws = 8 #experimental! (2021 Oct)
                        ws = 0.05 #dlambda (um)
                    else:
                        #ws = 3
                        ws = 0.1 #dlambda (um)
                    try:
                        #print('ve.processTemplate('+TPL_START+'.tpl,band="'"N"'"'+',gsmooth=0.15, wsmooth = '+'%d'%ws+',dophot=True,nbessel=0,resample=2,opdcutn=3)') #wsmmoth=1: no smoothing, 3: smooth without losing spectral resolution,  wsmooth=6 small loss of resolution
                        ###rdata = ve.processTemplate(TPL_START+'.tpl',band='N',gsmooth=0.15, wsmooth = ws,dophot=True,nbessel=0,resample=2,opdcutn=3) #old #wsmmoth=1: no smoothing, 3: smooth without losing spectral resolution,  wsmooth=6 small loss of resolution
                        
                        try:
                            rdata = ve.processTemplate(TPL_START+'.tpl',band='N',gsmooth=0.15,wsmooth=ws,dophot=False,opdcutn=3, switchframes=5) #new wsmooth = 0.025
                            ###rdata = wu.mrestore(RESDIR_BAND+'/'+TPL_START.replace(':', '_')+'.tpl.pk')
                            if(rdata is not None):
                                ofile = RESDIR_BAND+'/'+TPL_START.replace(':', '_')+'.tpl.pk'
                                wu.msave(ofile, rdata)
                            
                                #print(len(rdata),list(rdata[0].keys()))
                                #pk_2_oifits(RESDIR_BAND+'/'+TPL_START.replace(':', '_')+'.tpl.pk', RESDIR_BAND,oifits_template=PROGRAMDIR+'/template_oifits.fits',version='alpha') #old                            
                                print('Convert pk to oifits.')
                                pk_2_oifits(RESDIR_BAND+'/'+TPL_START.replace(':', '_')+'.tpl.pk', RESDIR_BAND,oifits_template=PROGRAMDIR+'/template_oifits.fits') #new
                                #print(rdata[0]['bcd1name'])
                                #fluxes = ve.templateAverageFlux(rdata, True)
                                #fluxes['header'] = rdata[0]['header']
                                #fluxes['bcd1name'] = ''
                                #fluxes['bcd2name'] = ''
                                #pk_2_oifits([fluxes], RESDIR_BAND,oifits_template=PROGRAMDIR+'/template_oifits.fits')
                        except FileNotFoundError as e:
                            print(e)
                    except IndexError as e:
                        print(e)
                os.chdir(current_dir)

    # ------------------------ make plots -----------------------
    if do_plot:
        for j in range(len(modes)):
            RESDIR = RESDIR_MAIN+'/'+modes[j]+'/' + NIGHT + '/' + TPL_START.replace(':', '_')
            if do_L:
                inputdir = RESDIR + '/Iter1/mat_raw_estimates.'+TPL_START.replace(':', '_')+'.HAWAII-2RG.rb/'
                # show_allred(inputdir,outputdir=inputdir+'/plots/',verbose=False,nbProc=6)
                show_allred_mosaic(inputdir, outputdir=inputdir + '/plots/',fn_pattern='RAW_INT',fit_model = False,annotate=False)
            if do_N:
                inputdir = RESDIR + '/Iter1/mat_raw_estimates.'+TPL_START.replace(':', '_')+'.AQUARIUS.rb/'
                # show_allred(inputdir, outputdir=inputdir + '/plots/',fn_pattern='', verbose=False,file_type='oifits',pro_catg='TARGET_RAW_INT',nbProc=6)
                # show_allred(inputdir, outputdir=inputdir + '/plots/',fn_pattern='', verbose=False,file_type='oifits',pro_catg='CALIB_RAW_INT',nbProc=6)
                if pipeline == 'ews':
                    show_allred_mosaic(inputdir, outputdir=inputdir + '/plots/',fn_pattern='RAW_INT',wl_lim=(7.3,13.4),fit_model = False,annotate=False)
                else:
                    show_allred_mosaic(inputdir, outputdir=inputdir + '/plots/',fn_pattern='RAW_INT',fit_model = False,annotate=False)

print('EXTERMINATE!')

if __name__ == "__main__":
    help(pk_2_oifits)
