from pathlib import Path
from astropy.io import fits as pyfits


def get_gravity_data(fits_file: Path)
    hdu = pyfits.open(file)

    # Read the data
    # index 10 = SC ; 20 = FT ; 11,12 = SC_POL ; 21,22 = FT_POL
    wave = hdu['OI_WAVELENGTH',10].data.field('EFF_WAVE')*1e6
    spectre = np.mean(hdu['OI_FLUX',10].data.field('FLUX'),0)
    visamp =  hdu['OI_VIS',10].data.field('VISAMP')
    visphi =  hdu['OI_VIS',10].data.field('VISPHI')
    closure = hdu['OI_T3',10].data.field('T3PHI')[:,:]
    #visampFT = hdu['OI_VIS',20].data.field('VISAMP')[:,3]
    ucoord = hdu['OI_VIS',10].data.field('UCOORD')
    vcoord = hdu['OI_VIS',10].data.field('VCOORD')

    # basename
    dicname = {i:n for i,n in zip(hdu['OI_ARRAY'].data.field('STA_INDEX'),
                                  hdu['OI_ARRAY'].data.field('STA_NAME'))}
    base = [dicname[i]+'-'+dicname[j]\
            for i, j in hdu['OI_VIS',10].data.field('STA_INDEX')]
    triplet = [dicname[i]+'-'+dicname[j]+'-'+dicname[k]\
               for i, j, k in hdu['OI_T3',20].data.field('STA_INDEX')]


if __name__ == "__main__":
    file = "C:/Users/perrautk/Documents/GRAVITY/GTO/HD100546/Myfile.fits"
