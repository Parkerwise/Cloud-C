'''
-*- coding: utf-8 -*-
Author: Parker Wise (Error and line counting by Ryan Cosgrove)
Date: 08-10-24
Description: quantifies richness towards each source
Python Version: 3.11.9
'''
from astropy import units as u
import pandas as pd
import pprint
import scipy
import numpy as np
from scipy import stats as st
from regions import Regions
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
def getError(freq, spectrum):
    binsize = 2
    bin_num = len(spectrum) / binsize
    binned_spectrum = st.binned_statistic(freq, spectrum, statistic='mean',
                                          bins=bin_num)
    sigma_bin = np.std(binned_spectrum.statistic)
    prelim_error = 4 * sigma_bin  # Four sigma of binned spectrum
    spectrum_ref = spectrum.copy()
    for channel in range(len(spectrum)):
        if spectrum_ref[channel].value > prelim_error or spectrum_ref[channel].value < -prelim_error:
            spectrum_ref[channel] = 0
    sigma_ref = np.std(spectrum_ref[0:-1].value)
    acceptable_error = 4 * sigma_ref  # Re-calculate error after first trimming
    spectrum_ref2 = spectrum_ref.copy()  # Copy spectrum and trim out values over error
    for channel in range(len(spectrum)):
        if spectrum_ref2[channel].value > acceptable_error or spectrum_ref2[channel].value < -acceptable_error:
            spectrum_ref2[channel] = 0
    sigma_ref2 = np.std(spectrum_ref2[0:-1].value)
    acceptable_error2 = 4 * sigma_ref2  # Re-calculate error after second trimming
    return acceptable_error2


def countLines(freq, spectrum, error):

    linechannels, properties = scipy.signal.find_peaks(spectrum, height=error,
                                                       prominence=error)
    linefreqs = freq[linechannels].value
    lineheights = spectrum[linechannels].value
    linenumber = len(linechannels)
    return linenumber, linefreqs, lineheights


def generateSpectra(cube, regionsList, regionIndex):
    subcube = cube.subcube_from_regions([regionsList[regionIndex]])
    spectrum = subcube.mean(axis=(1, 2))
    spectraError = getError(freq, spectrum)
    return spectrum, spectraError


# I decided to do dicts of dicts to keep track of multiple values
def linesTable(name, linefreqs, lineheights):
    dict = {
        f'{name} lines': {
            f'line {i:02}': {
                'central freq': freq,
                'line height': height
            } for i, (freq, height) in enumerate(zip(linefreqs, lineheights))
        }
    }
    return dict


'''
Part of our goal is to see how many lines identified in core c1 fall below
the noise of the spectra in other cubes. In order to do this we scale the lines
from c1 by the peak brightness of each source. These values are stored in
our catalog we generated in source-ID.py
'''
df = pd.read_csv('/home/pw/research/Cloud-C/results/tables/CloudC-catalog.csv')
peak_brightness = df.peak_brightness
brightness_ratios = [peak/peak_brightness[0] for peak in peak_brightness]


# Regions from which spectra is extracted
regions_file = "/home/pw/research/Cloud-C/fk5Regions.reg"
regpix = Regions.read(regions_file)
numberOfRegions = len(regpix)

c1Lines = {}
totalLines = {f"region {i}": 0 for i in range(numberOfRegions)}
totalabove = {f"region {i}": 0 for i in range(numberOfRegions)}
# region, num of lines
# Loop is run on each cube we're interested in
error = [0.042853665119683346, 0.11609788708836583, 0.034934569852619564,
         0.03281507998761665]
for i, (cube, name) in enumerate(zip(cubeList, cubeAbbreviation)):
    path = f"/home/pw/research/Cloud-C/co-data/{cube}"
    # # creates spectral cube object
    sc = SpectralCube.read(path)
    sc.allow_huge_operations = True

    # Converts units to GHz and K, feel free to change
    # #  defaults should be Hz and Jy/Beam
    sc = sc.with_spectral_unit(u.GHz)
    sc = sc.to(u.K)
    # # defines our frequency as a list
    freq, Dec, Ra = sc.world[:, 0, 0]

    # Calculate lines towards core C1
    c1Spectrum, c1Error = generateSpectra(sc, regpix, 0)
    c1NumOfLines, c1Freqs, c1Heights = countLines(freq, c1Spectrum, 4*error[i])
    c1LineProperties = linesTable(name, c1Freqs, c1Heights)
    c1Lines.update(c1LineProperties)
    totalLines["region 0"] += c1NumOfLines

    # C1 Heights in cube blank
    # scaled heights
    # scaled heights
    print(brightness_ratios)
    for j in range(1, numberOfRegions):
        scaledHeights = c1Heights * brightness_ratios[j]
        subcube = sc.subcube_from_regions([regpix[j]])
        spectrum = subcube.mean(axis=(1, 2))
        spectraError = error[i]*4
        numOfLines, linefreqs, linepeaks = countLines(freq, spectrum,
                                                      spectraError)
        freqRange = freq[-1].value-freq[0].value
        lineFrequency = numOfLines/freqRange
        totalLines[f"region {j}"] += numOfLines
        above = 0
        for height in scaledHeights:
            if height > spectraError:
                above += 1
        totalabove[f"region {j}"] += above
        print(name, f"region {j}:", f"above noise {above},",
              f'total lines {c1NumOfLines}')
totalabove["region 0"] = ""
print("\\begin{center}\n\\begin{tabular}{|c|c|c|}")
print("\\hline\nRegion & Lines Identified& Above Noise\\\\\n\\hlines")
regionNames = [8, 10, 12, 13, 16, 18, 20]
for i, region in enumerate(totalLines):
    print(f"{regionNames[i]} & {totalLines[region]} & {totalabove[region]} \\\\")
print("\\end{tabular}\n\\end{center}")
# scale lines, compare to noise
    # c1 Heights, 
# pprint.pprint(c1Lines)
# print(f'total number of lines in c1: {totalLines["region 00"]}')
# pprint.pprint(totalLines)
