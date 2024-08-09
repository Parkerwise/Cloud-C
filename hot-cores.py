'''
WIP
this code will make spectra for each ID'd source (see source-ID.py) and then
quantify it's line richness (lines/GHz) in the search for hot cores which are
line rich
'''
# code yanked from spectra script
from astropy import units as u
from regions import Regions
import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
import sys
import regions
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# Reading in files
# # This should be the file where we're extracting spectra
path = "/home/pw/research/Cloud-C/co-data/B.Dust_Ridge_sci.spw31.cube.I.pbcor.fits"
# # creates spectral cube object
sc = SpectralCube.read(path)
sc.allow_huge_operations = True


# Converts units to GHz and K, feel free to change
# #  defaults should be Hz and Jy/Beam
sc = sc.with_spectral_unit(u.GHz)
sc = sc.to(u.K)
# # defines our frequency as a list
freq, Dec, Ra = sc.world[:, 0, 0]


# defines a subcube around your source
# # From regions file

regions_file = "/home/pw/research/Cloud-C/co-data/c1-by-eye"
regpix = Regions.read(regions_file)
subcube = sc.subcube_from_regions([regpix[0]])

# averages values within aperture for each channel
spectrum = subcube.mean(axis=(1, 2))


fig1 = plt.figure(1, figsize=(15, 2), dpi=250)
plt.plot(freq, spectrum, lw=1, drawstyle='steps-mid', color="SteelBlue",
         label="spectrum")

# Calculating and plotting noise, specific to your data
# # Noise_upper to Noise_lower should be line free
Noise_upper = 825  # defines some line free channels to calculate STD
Noise_lower = 800
sigma = np.std(spectrum[Noise_lower:Noise_upper].value)
three_sigma = 3*sigma
plt.hlines(three_sigma, freq[0].value, freq[1916].value, colors="red",
           label='', ls="--")
plt.hlines(-three_sigma, freq[0].value, freq[1916].value, colors="red",
           label='', ls="--")
plt.xlim(freq[0].value, freq[1916].value)
plt.fill_between(freq.value, three_sigma, -three_sigma, alpha=0.2,
                 color='red', label='Error')

# Formatting plot
fig1.supxlabel("Frequency (GHz)", fontsize=10)
fig1.supylabel('Brightness Temp. (K)', fontsize=10)
plt.legend(fontsize=14, loc="upper left")
plt.tight_layout()  # Adjust params to avoid overlap and decrease white space
# good to save as both png and pdf
'''
plt.savefig("spectra.png")
plt.savefig("spectra.pdf")
'''
plt.show()
# Saving figures

