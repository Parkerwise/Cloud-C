'''
-*- coding: utf-8 -*-
Author: Parker Wise
Date:
Description: make spectra for ID'd sources (see source-ID.py)
Python Version: 3.11.9
'''
# code yanked from spectra script
from astropy import units as u
from regions import Regions
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
    '''
    # Noise is not included at this moment, channels are cube dependent
    Noise_upper = 825  # defines some line free channels to calculate STD
    Noise_lower = 800
    '''

    fig1 = plt.figure(1, figsize=(15, 2*numberOfRegions), dpi=250)
    for i in range(numberOfRegions):
        subcube = sc.subcube_from_regions([regpix[i]])
        spectrum = subcube.mean(axis=(1, 2))
        ax1 = plt.subplot(numberOfRegions, 1, i+1)
        ax1.plot(freq, spectrum, lw=1, drawstyle='steps-mid', color="SteelBlue")
        ax1.set_title(f"region {i}")
        '''
        sigma = np.std(spectrum[Noise_lower:Noise_upper].value)
        three_sigma = 3*sigma
        plt.hlines(three_sigma, freq[0].value, freq[1916].value, colors="red",
                   label='', ls="--")
        plt.hlines(-three_sigma, freq[0].value, freq[1916].value, colors="red",
                   label='', ls="--")
        plt.xlim(freq[0].value, freq[1916].value)
        plt.fill_between(freq.value, three_sigma, -three_sigma, alpha=0.2,
                         color='red', label='Error')
        '''
    fig1.supxlabel("Frequency (GHz)", fontsize=10)
    fig1.supylabel('Brightness Temp. (K)', fontsize=10)
    plt.tight_layout()  # Adjust params to avoid overlap and decrease white space
    # good to save as both png and pdf
    plt.savefig(f"/home/pw/research/Cloud-C/results/spectra/{name}/{name}-line-richness.png")
    plt.savefig(f"/home/pw/research/Cloud-C/results/spectra/{name}/{name}-line-richness.pdf")
