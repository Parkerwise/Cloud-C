'''
-*- coding: utf-8 -*-
Author: Parker Wise (Error and line counting by Ryan Cosgrove)
Date: 08-10-24
Description: quantifies richness towards each source
Python Version: 3.11.9
'''
from astropy import units as u
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


def linesTable(name, linefreqs, lineheights):
    dict = {
        f'{name} lines': {
            f'line {i}': {
                'central freq': freq,
                'line height': height
            } for i, (freq, height) in enumerate(zip(linefreqs, lineheights))
        }
    }
    return dict


c1Lines = {}
# Loop is run on each cube we're interested in
for cube, name in zip(cubeList, cubeAbbreviation):
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

    # Regions from which spectra is extracted
    regions_file = "/home/pw/research/Cloud-C/fk5Regions.reg"
    regpix = Regions.read(regions_file)
    numberOfRegions = len(regpix)

    # Calculate lines towards core C1
    c1Spectrum, c1Error = generateSpectra(sc, regpix, 1)
    c1NumOfLines, c1Freqs, c1Heights = countLines(freq, c1Spectrum, c1Error)
    c1LineProperties = linesTable(name, c1Freqs, c1Heights)
    c1Lines.update(c1LineProperties)

    for i in range(numberOfRegions):
        subcube = sc.subcube_from_regions([regpix[i]])
        spectrum = subcube.mean(axis=(1, 2))
        spectraError = getError(freq, spectrum)
        numOfLines, linefreqs, linepeaks = countLines(freq, spectrum,
                                                      spectraError)
        freqRange = freq[-1].value-freq[0].value
        lineFrequency = numOfLines/freqRange
