'''
-*- coding: utf-8 -*-
Author: Parker Wise
Date: 08-10-24
Description: make spectra for ID'd sources (see source-ID.py)
Python Version: 3.11.9
'''
# code yanked from spectra script
from astropy import units as u
import numpy as np
from scipy import stats as st
from regions import Regions
import matplotlib.pyplot as plt
import pyspeckit as psk
from spectral_cube import SpectralCube
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# Reading in files
# # This should be the file where we're extracting spectra
cube_list = [
    "A.Dust_Ridge_12C18O.cube.I.pbcor.fits",
    "A.Dust_Ridge_13C16O.cube.I.pbcor.fits",
    "B.Dust_Ridge_12C16O_1-0.cube.I.pbcor.fits",
    "B.Dust_Ridge_12C17O.cube.I.pbcor.fits",
]

cube_abbreviation = [
    "A.SPW25",  # 12C18O
    "A.SPW27",  # 13C16O
    "B.SPW25",  # 12C16O
    "B.SPW27",  # 12C17O
]

isotopologue_range = [
    (109.88, 109.90),  # 12C18O
    (110.15, 110.25),  # 13C16O
    (115.20, 115.35),  # 12C16O
    (112.30, 112.40),  # 12C17O
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
    # Copy spectrum and trim out values over error
    spectrum_ref2 = spectrum_ref.copy()
    for channel in range(len(spectrum)):
        if spectrum_ref2[channel].value > acceptable_error or spectrum_ref2[channel].value < -acceptable_error:
            spectrum_ref2[channel] = 0
    sigma_ref2 = np.std(spectrum_ref2[0:-1].value)
    acceptable_error2 = 4 * sigma_ref2  # Re-calculate error after second trimming
    return acceptable_error2


error = [
    0.042853665119683346,
    0.11609788708836583,
    0.034934569852619564,
    0.03281507998761665
]
# Loop is run on each cube we're interested in
for cube, name, error, bounds in zip(cube_list, cube_abbreviation,
                                     error, isotopologue_range):
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

    # fig, axs = plt.subplots(numberOfRegions/2, 2, sharex=True,
    #                         figsize=(6, numberOfRegions), dpi=250,
    #                         layout='constrained')
    names = [1, 7, 2, 4, 5, 3, 6]
    fig1 = plt.figure(1, figsize=(10, 10), dpi=250)
    for i in range(numberOfRegions):
        subcube = sc.subcube_from_regions([regpix[i]])
        spectrum = subcube.mean(axis=(1, 2))
        order = names[i]
        '''
        print(order)
        # axs[order].tick_params(axis='both', which='major', labelsize=10)
        axs[order].plot(freq, spectrum, lw=1,
                        drawstyle='steps-mid', color="SteelBlue")
        axs[order].set_title(f"Region {names[i]}")
        axs[order].set_xlim(freq[0].value, freq[-1].value)
        axs[order].set_ylim(-0.5, 0.5)
        spectraError = error*4
        axs[order].hlines(spectraError, freq[0].value, freq[-1].value,
                          colors="red", ls="--")
        axs[order].hlines(-spectraError, freq[0].value, freq[-1].value,
                          colors="red", ls="--")
        axs[order].fill_between(freq.value, spectraError, -spectraError, alpha=0.2,
                                color='red', label='Error')
        '''
        zoom = np.where((freq.value > bounds[0])*(freq.value < bounds[1]))
        zoom_12C180 = spectrum[zoom]
        zoom_freq = freq[zoom]
        zeros = np.zeros(np.size(zoom))
        wheremeas = np.where((zoom_freq.value > 115.33)*(zoom_freq.value < 115.35))
        error = zeros + error
        sp = psk.Spectrum(
            data=zoom_12C180, xarr=zoom_freq, error=error, unit='K'
        )
        sp.plotter(axis=plt.subplot(4, 2, order),
                   label=f"Reg. {order}",
                   ymin=-1, ymax=3
                   )
        sp.plotter.label(xlabel="", ylabel="")
        # sp.specfit.plot_fit(annotate=False)
        sp.specfit(fittype='gaussian', annotate=False)
        # sp.plotter.savefig('12C16Ogaussianprofileplot.pdf')
        # plt.xticks(rotation=45)
        plt.legend()
    print(sp.specfit.parinfo)
    fig1.suptitle(name)
    # plt.title(name, fontsize=10)
    # fig.subplots_adjust(wspace=0.1,
    #                     hspace=0.3)
    # fig.supylabel('Brightness Temp. (K)', fontsize=10)
    # # good to save as both png and pdf
    #
    fig1.subplots_adjust(wspace=0, hspace=0)
    fig1.tight_layout()
    plt.savefig(f"foo-{name}.png", dpi=250)
    plt.clf()  # Removes this cubes spectra so it won't be plotted on the next
