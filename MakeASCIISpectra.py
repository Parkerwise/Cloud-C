'''
-*- coding: utf-8 -*-
Author: Parker Wise
Date: 2024-07-05
Description: outputs Spectra as ASCII readable by xclass
Python Version: 3.11.9
'''

from astropy import units as u
from spectral_cube import SpectralCube
import sys
import regions
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# This should be the file where we're extracting spectra
path = input("FITS File: \n")

# the center of the source in pixels from where we are extracting spectra
central_xcoord = int(input("Central X Coord (px): \n"))
central_ycoord = int(input("Central y coord (px): \n"))

# radius of the source in pixels
radius = int(input("Radius (px): \n"))
csv_name = input("Name of new spectra file name: \n")
subgroups = int(input("Split your spectra into how many files?: \n"))
# reads in file
sc = SpectralCube.read(path)
sc.allow_huge_operations = True

# Converts units to MHz and K, feel free to change
# defaults should be Hz and Jy/Beam
sc = sc.with_spectral_unit(u.MHz)
sc = sc.to(u.K)

# defines our frequency as a variable
freq, Dec, Ra = sc.world[:, 0, 0]

# defines a subcube around your source
regpix = regions.CirclePixelRegion(regions.PixCoord(central_xcoord,
                                                    central_ycoord), radius=radius)
subcube = sc.subcube_from_regions([regpix])

# averages values within subcube
spectrum = subcube.mean(axis=(1, 2))

for i in range(subgroups):
    with open(f"{csv_name}.MHz.{i}.dat", "w") as spectratext:
        for j in range(round(len(freq)/9*i), round(len(freq)/9*(i+1))):
            spectratext.write(f"\t{freq[j].value}\t{spectrum[j].value}\n")
with open(f"{csv_name}.MHz.complete.dat", "w") as spectratext:
    for j in range(len(freq)):
        spectratext.write(f"\t{freq[j].value}\t{spectrum[j].value}\n")
