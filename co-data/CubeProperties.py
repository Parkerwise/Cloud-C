#work in progress
import astropy.io.fits as fits #6.1.0
import numpy as np
from spectral_cube import SpectralCube
#from astropy.coordinates import SkyCoord
#from astropy import units as u  

#plotting and stuff
#formatting stuff
path ="A.Dust_Ridge_12C18O.cube.I.pbcor.fits"
image=fits.getdata(path)
def cube_properties(path):
    header=fits.getheader(path)
    name=path[0]+'.'+header['FILNAM04']
    xdim=header['CTYPE1']
    ydim=header['CTYPE2']
    pixelScale=round(float(header['CDELT2'])*3600,3)
    channels=len(image[0])
    #freqRange=header['']
    #freqWidth=header['']
    #bmaj=header['']
    #bmin=header['']
    #bpa=header['']
    #spatialResolutions=header['']
    #spectralResolutions=header['']
    #channelWidth=header['']
    return f"{name},{xdim}, {ydim}, {pixelScale}, {channels}"#, {freqRange}, {freqWidth}, {bmaj}, {bmin}, {bpa},{spatialResolutions}, {spectralResolutions}, {channelWidth}"
with open("cube-properties.csv",'w') as table:
    table.write(cube_properties(path))
