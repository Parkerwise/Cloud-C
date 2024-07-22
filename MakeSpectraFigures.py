'''
-*- coding: utf-8 -*-
Author: Parker Wise
Date: 2024-07-01
Description: Makes spectra with set window sizes
Python Version: 3.11.9
'''
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy import units as u
import pylab
import math
from regions import Regions
import numpy as np
from spectral_cube import SpectralCube
import sys
from matplotlib.ticker import (MultipleLocator)
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

path = input("FITS Path:\n")
width_per_plot = input("Desired width per plot (in GHz): \n")
regions_file = input("Regions File Path: \n")
filename = input("filename:\n")
header = fits.getheader(path)
w1 = WCS(header)

sc = SpectralCube.read(path)
sc.allow_huge_operations = True
sc_Ghz = sc.with_spectral_unit(u.GHz)
sc_Ghz = sc_Ghz.to(u.K)
freq, Dec, Ra = sc_Ghz.world[:, 0, 0]

regpix = Regions.read(regions_file)
subcube = sc_Ghz.subcube_from_regions([regpix[0]])
spectrum = subcube.mean(axis=(1, 2))

channel_width = freq[1].value-freq[0].value
full_width = freq[-1].value-freq[0].value
channels_per_window = math.ceil(width_per_plot/channel_width)
number_of_windows = math.ceil(full_width/width_per_plot)

fig1 = pylab.figure(1, figsize=(15, 2*number_of_windows))
for i in range(number_of_windows):
    if i+1 != number_of_windows:
        ax1 = pylab.subplot(number_of_windows, 1, i+1)
        ax1.plot(freq[channels_per_window*i:channels_per_window*(i+1)],
                 spectrum[channels_per_window*i:channels_per_window*(i+1)],
                 lw=1, drawstyle='steps-mid', color="SteelBlue")
        pylab.xlim(freq[channels_per_window*i].value,
                   freq[channels_per_window*(i+1)].value)
        ax1.xaxis.set_major_locator(MultipleLocator(.05))
        ax1.xaxis.set_minor_locator(MultipleLocator(.01))
    else:
        ax1 = pylab.subplot(number_of_windows, 1, number_of_windows)
        empty_size = channels_per_window*(i+1)-len(freq)
        empty_list = [0]*empty_size
        extra_frequencies = [(j+1)*channel_width+freq[-1].value for j in range(empty_size)]
        remaining_freqs = freq[channels_per_window*i:].value
        freq_with_zeros = np.append(remaining_freqs, extra_frequencies)
        remaining_data = spectrum[channels_per_window*i:]
        data_with_zeros = np.append(remaining_data, empty_list)
        ax1.plot(freq_with_zeros, data_with_zeros.value, lw=1,
                 drawstyle='steps-mid', color="SteelBlue")
        pylab.xlim(freq_with_zeros[0], freq_with_zeros[-1])
        ax1.xaxis.set_major_locator(MultipleLocator(.05))
        ax1.xaxis.set_minor_locator(MultipleLocator(.01))
fig1.supxlabel("Frequency (GHz)", fontsize=10)
fig1.supylabel('Brightness Temp. (K)', fontsize=10)
plt.rcParams['text.usetex'] = True
fig1.tight_layout(rect=[0.005, 0, 1, 0.001])
'''
plt.savefig(f"{filename}.divided.pdf")
plt.savefig(f"{filename}.divided.png")
'''
