'''
-*- coding: utf-8 -*-
Author: Parker Wise
Date: 07-21-2024
Description: ID's structures in Cloud c
Python Version: 3.11.9
'''
# librarys
from radio_beam import Beam  # 0.3.7
from astrodendro.plot import DendrogramPlotter
from astrodendro import Dendrogram, pp_catalog
import astropy.io.fits as fits  # 6.1.0
# astrodendro viewer not compatible with matplotlib >= 3.7
import matplotlib.pyplot as plt  # 3.6.0
from astropy.wcs import WCS
from astropy import units as u
import numpy as np  # 1.26.4
import reproject  # 0.13.1
from reproject.mosaicking import find_optimal_celestial_wcs
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# Importing image
path = "/home/pw/research/Cloud-C/co-data/CloudC_mJy.fits"
image = fits.getdata(path)
header = fits.getheader(path)
w1 = WCS(header)
w1 = w1.dropaxis(3)
w1 = w1.dropaxis(2)
image_2D = np.squeeze(image)


# masking without using regions library
def addMask(image, radius, center, remove_inside=True):
    Y, X = np.ogrid[:350, :350]
    dist = np.sqrt((X-center[0])**2+(Y-center[1])**2)
    mask = dist <= radius
    if remove_inside is True:
        for (i, j), bool in np.ndenumerate(mask):
            if bool:
                image[i, j] = np.nan
    if remove_inside is False:
        for (i, j), bool in np.ndenumerate(mask):
            if not bool:
                image[i, j] = np.nan


addMask(image_2D, 110, (166, 175), remove_inside=False)  # crops noisy edges

# converts image to galactic coordinates
wcs_out, shape_out = find_optimal_celestial_wcs([(image_2D, w1)],
                                                frame='galactic')
imageGalactic, c_footprint = reproject.reproject_interp((image_2D, w1),
                                                        wcs_out,
                                                        shape_out=shape_out)

# Computing Dendrogram
sigma = 0.41  # mJy/beam, std of noise
cloudDendrogram = Dendrogram.compute(imageGalactic, min_value=sigma,
                                     min_delta=2*sigma, wcs=wcs_out)

# Creates catalog of dendrogram structures
catalogMetadata = {}
catalogMetadata['data_unit'] = u.mJy / u.beam
catalogMetadata['spatial_scale'] = 0.28 * u.arcsec
catalogMetadata['beam_major'] = 2.259 * u.arcsec  # FWHM
catalogMetadata['beam_minor'] = 1.590 * u.arcsec  # FWHM
cloudCatalog = pp_catalog(cloudDendrogram, catalogMetadata)
cloudCatalog.write('/home/pw/research/Cloud-C/results/tables/CloudC-catalog.csv',
                   format='ascii.csv', overwrite=True)
cloudCatalog.write('/home/pw/research/Cloud-C/results/tables/CloudC-catalog.ecsv',
                   format='ascii.ecsv', overwrite=True, delimiter=',')
'''
catalog can be saved many different ways. See ASCII tables in astropy docs
cat.write('CloudC-catalog.fits', format='fits', overwrite=True)
cat.write('CloudC-catalog.tex', format='ascii.latex', overwrite=True)
'''
cloudPlots = DendrogramPlotter(cloudDendrogram)  # plotter instance

# Dendrogram Plot
dendrogramFigure = plt.figure(0, figsize=(14, 14), constrained_layout=True)
dendrogramAx = dendrogramFigure.add_subplot(111)
cloudPlots.plot_tree(dendrogramAx)
dendrogramAx.set_xlabel("Structure")
dendrogramAx.set_ylabel("Flux")
dendrogramAx.set_ylim(0, 15)
plt.savefig("/home/pw/research/Cloud-C/results/continuum/CloudC-dendrogram.pdf")
plt.savefig("/home/pw/research/Cloud-C/results/continuum/CloudC-dendrogram.png")

# continuum and contour plot
continuumFigure = plt.figure(1, figsize=(14, 14), constrained_layout=True)
continuumAx = continuumFigure.add_subplot(111, projection=wcs_out)
'''
structures were plotted as contours
each dendrogram level is associated with a color
if dendrogram has more than 6 levels, more colors must be specified
this same affect can be achieved with many levels if you use a color map
this was done previously in commit 564d8e2
'''
colors = ["red", "orange", "yellow", "green", "cyan", "purple"]
for structure in cloudDendrogram:
    currentLevel = structure.level
    cloudPlots.plot_contour(continuumAx, structure=structure,
                            color=colors[currentLevel])
lon = continuumAx.coords[0]
lon.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
lat = continuumAx.coords[1]
lat.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
plt.xlim(80, 410)
plt.ylim(80, 410)

# plots beam
my_beam = Beam.from_fits_header(header)
ycen_pix, xcen_pix = 100, 375  # location of beam on plot
pixscale = 0.28 * u.arcsec
ellipse_artist = my_beam.ellipse_to_plot(xcen_pix, ycen_pix, pixscale)
continuumImage = plt.imshow(imageGalactic, cmap='Greys_r', vmax=5)
plt.gca().add_patch(ellipse_artist)  # plots beam
ellipse_artist.set_facecolor("white")
ellipse_artist.set_edgecolor("black")


# plots scale bar
x = [90, 180]
y = [90, 90]
plt.plot(np.array([x[0], x[1]]), np.array([y[0], y[1]]),
         color="black", linewidth=3)
# scale is calculated using small angle approximation
# to get angular seperation between pixels, use this function
'''
def angularSeperation(beginningPixel, endingPixel):
    beginningPixel = w1.pixel_to_world(x[0], y[0])
    endingPixel = w1.pixel_to_world(x[1], y[1])
    return beginningPixel.separation(endingPixel)
'''
plt.text(125, 100, '1pc', fontsize=14, color='black')

# formats plot
lon.set_ticks(size=-3)
lat.set_ticks(size=-3)
plt.xlabel('Galactic Longtitude', fontsize=20, labelpad=1)
plt.ylabel('Galactic Latitutude', fontsize=20, labelpad=1)
continuumAx.tick_params(axis='both', which='major', labelsize=15)
plt.annotate('Continuum', fontsize=15, xy=(0.02, 0.91), xycoords="axes fraction")
plt.savefig("/home/pw/research/Cloud-C/results/continuum/CloudC-structure-contours.pdf")
plt.savefig("/home/pw/research/Cloud-C/results/continuum/CloudC-structure-contours.png")