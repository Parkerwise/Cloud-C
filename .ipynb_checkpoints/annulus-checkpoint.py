import astropy.io.fits as fits
import matplotlib.pyplot as plt    
from astropy.wcs import WCS                 
from astropy import units as u  
import pylab
import regions
import numpy as np
from spectral_cube import SpectralCube
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
fig, (ax1, ax2)=plt.subplots(1,2)
#Takes pixel position in Cont --> pixel position in new cube
path="B.Dust_Ridge_sci.spw25_27_29_31.cont.I.tt0.pbcor.fits"
header=fits.getheader(path)
w1=WCS(header)
w1 = w1.dropaxis(3)
w1 = w1.dropaxis(2)
#finds the pixel positions for the RA and dec
#find RA DEC of x,y,z
core_coords=w1.pixel_to_world(175,166)
#reads in file
path="B.Dust_Ridge_12C16O_1-0.cube.I.pbcor.fits"
#find the corresponding pixel positions for those coords in new image
header=fits.getheader(path)
w1=WCS(header)
image_coords=core_coords.to_pixel(w1,0,mode="wcs")
#makes figure
fig1=pylab.figure(1,figsize=(15,2),dpi=250)
#plots spectrum
sc=SpectralCube.read(path)
sc.allow_huge_operations=True 
sc_Ghz=sc.with_spectral_unit(u.GHz)
sc_Ghz=sc_Ghz.to(u.K)
freq,Dec,Ra = sc_Ghz.world[:,0,0]
print(image_coords)


ogsubcube=sc_Ghz.hdu.data[:,170:180,161:171]
ogspectrum = np.mean(ogsubcube,axis=(1,2))


regpix = regions.RectanglePixelRegion(regions.PixCoord(166,175), width=10, height=10)  
subcube = sc_Ghz.subcube_from_regions([regpix])  
spectrum = subcube.mean(axis=(1, 2))
diff = [ogspectrum[i]-spectrum[i].value for i in range(len(spectrum))]
print(f"ogsubcube:\n {np.shape(ogsubcube)}")
print(f"subcube:\n {subcube}")


ax1.plot(freq,diff,lw=1,drawstyle='steps-mid',color="SteelBlue",label="diff")
ax2.plot(freq,ogspectrum,lw=1,drawstyle='steps-mid',color="SteelBlue",label="OGSpectrum")
ax2.plot(freq,spectrum,lw=1,drawstyle='steps-mid',color="green",label="spectrum")
pylab.xlabel("Frequency (GHz)", fontsize=10) 
pylab.ylabel('Brightness Temp. (K)',fontsize=10)
plt.rcParams['text.usetex'] = True
plt.legend(fontsize=14,loc="upper left")
plt.show()
