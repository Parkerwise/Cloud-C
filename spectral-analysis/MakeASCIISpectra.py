#for questions please reach out to wiseparker03@ku.edu
#this code will produce a CSV of your spectra (in MHz) given a data cube and a source
import astropy.io.fits as fits
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from astropy.wcs import WCS                 
from astropy import units as u  
import numpy as np
from spectral_cube import SpectralCube
import sys
import regions
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

#This should be the file where we're extracting spectra
path=input("FITS File: \n")

#the center of the source in pixels from where we are extracting spectra
central_xcoord=int(input("Central X Coord (px): \n"))
central_ycoord=int(input("Central y coord (px): \n"))

#radius of the source in pixels
radius=int(input("Radius (px): \n"))


csv_name=input("Name of new spectra file name: \n")

#reads in file
sc=SpectralCube.read(path)
sc.allow_huge_operations=True 

#Converts units to MHz and K, feel free to change
#defaults should be Hz and Jy/Beam
sc=sc.with_spectral_unit(u.MHz)
sc=sc.to(u.K)

#defines our frequency as a variable
freq,Dec,Ra = sc.world[:,0,0]

#defines a subcube around your source
#166,175, 8
regpix = regions.CirclePixelRegion(regions.PixCoord(central_xcoord,central_ycoord), radius=radius)  
subcube = sc.subcube_from_regions([regpix])  

#averages values within subcube
spectrum = subcube.mean(axis=(1, 2))


with open(f"{csv_name}", "w") as spectratext:
    for i in range(len(freq)):
        spectratext.write(f"\t{freq[i].value}\t{spectrum[i].value}\n")
plt.plot(freq,spectrum)
plt.show()


