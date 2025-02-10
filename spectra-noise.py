'''
-*- coding: utf-8 -*-
Author: Parker Wise
Date: 08-10-24
Description: make spectra for ID'd sources (see source-ID.py)
Python Version: 3.11.9
'''
# code yanked from spectra script
import astropy.io.fits as fits  # 6.1.0
from astropy import units as u
from astropy.coordinates import Angle
from astropy.wcs import WCS
import numpy as np
import regions
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# Reading in files
# # This should be the file where we're extracting spectra
cubeList = [
    'A.Dust_Ridge_sci.spw29.cube.I.pbcor.fits',  # A29
    'A.Dust_Ridge_sci.spw31.cube.I.pbcor.fits',  # A31
    'B.Dust_Ridge_sci.spw29.cube.I.pbcor.fits',  # B29
    'B.Dust_Ridge_sci.spw31.cube.I.pbcor.fits',  # B31
]
# For naming final products/saving into results folder
cubeAbbreviation = [
    "A.spw29",  # A29
    "A.spw31",  # A31
    "B.spw29",  # B29
    "B.spw31",  # B31
]


# Define a function to obtain the error from a given spectrum

# Loop is run on each cube we're interested in
noiseLimits = [
              [1050, 1900],  # A29
              [800, 1150],  # A31
              [100, 600],  # B29
              [1300, 1900]]  # B31
fig1 = plt.figure(1, figsize=(30, 10))
for i, (cube, name) in enumerate(zip(cubeList, cubeAbbreviation)):
    path = f"/home/pw/research/Cloud-C/co-data/{cube}"
    header = fits.getheader(path)
    w1 = WCS(header)
    w1 = w1.dropaxis(3)
    w1= w1.dropaxis(2)
    # # creates spectral cube object
    sc = SpectralCube.read(path)
    sc.allow_huge_operations = True

    # Converts units to GHz and K, feel free to change
    # #  defaults should be Hz and Jy/Beam
    sc = sc.with_spectral_unit(u.GHz)
    sc = sc.to(u.K)
    # # defines our frequency as a list
    freq, Dec, Ra = sc.world[:, 0, 0]

    # Regions from which spectra is extracted


    #######
    bmaj = header['BMAJ'] * 3600 /0.28
    bmin = header['BMIN'] * 3600 /0.28
    bpa = header['BPA']
    central_xcoord = 182
    central_ycoord = 195
    pixCoord = regions.PixCoord(central_xcoord, central_ycoord)
    skyCoord = w1.pixel_to_world(central_xcoord, central_ycoord)
    regpix = regions.EllipsePixelRegion(center=pixCoord, width=bmaj,
                                        height=bmin, angle=Angle(bpa, 'deg'))

    ######
    subcube = sc.subcube_from_regions([regpix])
    spectrum = subcube.mean(axis=(1, 2))
    ax1 = plt.subplot(4, 1, i+1)
    ax1.plot(freq, spectrum, lw=1, drawstyle='steps-mid', color="SteelBlue")
    ax1.set_xlim(freq[0].value, freq[-1].value)
    ax1.set_ylim(-0.5, 0.5)
    noiseUpper = noiseLimits[i][1]
    noiseLower = noiseLimits[i][0]
    plt.vlines(freq[noiseUpper].value, -1, 2, colors="black")
    plt.vlines(freq[noiseLower].value, -1, 2, colors="black")
    sigma = np.std(spectrum[noiseLower:noiseUpper].value)
    print(name, sigma)
    plt.hlines(sigma, freq[0].value, freq[-1].value,
               colors="red", ls="--")
    plt.hlines(-sigma, freq[0].value, freq[-1].value,
               colors="red", ls="--")
    plt.fill_between(freq.value, sigma, -sigma, alpha=0.2,
                     color='red', label='Error')
    fig1.supxlabel("Frequency (GHz)", fontsize=10)
    fig1.supylabel('Brightness Temp. (K)', fontsize=10)
    plt.tight_layout()  # Adjust params to avoid overlap and decrease white space
    # good to save as both png and pdf
    # plt.savefig(f"/home/pw/research/Cloud-C/results/spectra/{name}/{name}-line-richness.png")
    # plt.savefig(f"/home/pw/research/Cloud-C/results/spectra/{name}/{name}-line-richness.pdf")
plt.show()
