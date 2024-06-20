#foo
import astropy.io.fits as fits
import matplotlib.pyplot as plt    
from astropy.wcs import WCS                 
from astropy import units as u  
import pylab
import math
import pyspeckit
from regions import Regions
import numpy as np
from spectral_cube import SpectralCube
import sys
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

path="B.Dust_Ridge_sci.spw31.cube.I.pbcor.fits"
regions_file="c1-by-eye"
header=fits.getheader(path)
w1=WCS(header)

sc=SpectralCube.read(path)
sc.allow_huge_operations=True 
sc_Ghz=sc.with_spectral_unit(u.GHz)
sc_Ghz=sc_Ghz.to(u.K)
freq,Dec,Ra = sc_Ghz.world[:,0,0]

regpix = Regions.read(regions_file)
subcube = sc_Ghz.subcube_from_regions([regpix[0]])  
spectrum = subcube.mean(axis=(1, 2))


Noise_upper=825
Noise_lower=800
sigma=np.std(spectrum[Noise_lower:Noise_upper].value)
three_sigma=3*sigma
zoom = np.where((freq.value>102.4)*(freq.value<102.6))
zoom_spectrum=spectrum[zoom]
zoom_freq = freq[zoom]
#Calculates the error 
error = np.zeros(np.size(zoom))
meas = sigma #uses sigma calc from previous block
error = error + meas

#plots the gaussian fit
sp = pyspeckit.Spectrum(data=zoom_spectrum,xarr=zoom_freq,error=error,unit='K')
fig1 = pylab.figure(1,figsize=(5,5),dpi=250)
sp.plotter(axis=pylab.subplot(1,1,1),title="B.SPW31")
sp.specfit(fittype='gaussian')
#sp.plotter.savefig('B.spw31.brightest.pdf')
#sp.plotter.savefig('B.spw31.brightest.png')
plt.show()
print(sp.specfit.parinfo)
