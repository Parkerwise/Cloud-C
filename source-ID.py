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
from astrodendro.analysis import PPStatistic
from astrodendro import Dendrogram, pp_catalog
import astropy.io.fits as fits  # 6.1.0
from regions import EllipseSkyRegion
# astrodendro viewer not compatible with matplotlib >= 3.7
import matplotlib.pyplot as plt  # 3.6.0
from astropy.wcs import WCS
from astropy import units as u
import numpy as np  # 1.26.4
import reproject  # 0.13.1
import matplotlib
from reproject.mosaicking import find_optimal_celestial_wcs
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# Importing image
matplotlib.rcParams.update({'font.size': 10})
plt.rcParams['xtick.major.pad'] = '8'
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


addMask(image_2D, 100, (165, 190), remove_inside=False)  # crops noisy edges

# converts image to galactic coordinates
wcs_out, shape_out = find_optimal_celestial_wcs([(image_2D, w1)],
                                                frame='galactic')
imageGalactic, c_footprint = reproject.reproject_interp((image_2D, w1),
                                                        wcs_out,
                                                        shape_out=shape_out)


# Computing Dendrogram
sigma = 0.36  # mJy/beam, std of noise
cloudDendrogram = Dendrogram.compute(imageGalactic, min_value=sigma,
                                     min_delta=1*sigma, wcs=wcs_out)
# manually removes unwanted trunk and leaves
unwantedIDS = [22, 24, 14, 7, 23, 21, 19, 17]
unwantedStructs = [24, 0, 1, 2, 5, 22, 21]
unwantedLeaves = [13, 8, 6, 0, 12, 10, 11]
'''
cloudDendrogram.trunk.remove(cloudDendrogram.trunk[0])
cloudDendrogram.trunk.remove(cloudDendrogram.trunk[1])
cloudDendrogram.trunk[0].descendants.remove(cloudDendrogram.leaves[12])
cloudDendrogram.trunk[0].descendants.remove(cloudDendrogram.leaves[10])
cloudDendrogram.trunk[0].descendants.remove(cloudDendrogram.leaves[0])
print(cloudDendrogram.trunk[0].descendants)
'''
# Creates catalog of dendrogram structures
catalogMetadata = {}
catalogMetadata['data_unit'] = u.mJy / u.beam
catalogMetadata['spatial_scale'] = 0.28 * u.arcsec
catalogMetadata['beam_major'] = 2.259 * u.arcsec  # FWHM
catalogMetadata['beam_minor'] = 1.590 * u.arcsec  # FWHM
catalogMetadata['wcs'] = wcs_out
cloudCatalog = pp_catalog(cloudDendrogram.leaves, catalogMetadata)
cloudCatalogAdjusted = cloudCatalog
for struct in unwantedIDS:
    catalogMask = (cloudCatalogAdjusted['_idx'] != struct)
    cloudCatalogAdjusted = cloudCatalogAdjusted[catalogMask]
cloudCatalogAdjusted['peak_brightness'] = 0. * u.mJy / u.beam
j = 0
for i, leaf in enumerate(cloudDendrogram.leaves):
    if i not in unwantedLeaves:
        cloudCatalogAdjusted['peak_brightness'][j] = leaf.get_peak()[1]
        j += 1
print(cloudCatalogAdjusted)
cloudCatalogAdjusted.write('/home/pw/research/Cloud-C/results/tables/CloudC-catalog.csv',
                   format='ascii.csv', overwrite=True)
cloudCatalogAdjusted.write('/home/pw/research/Cloud-C/results/tables/CloudC-catalog.ecsv',
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
contourFigure = plt.figure(1, figsize=(14, 14), constrained_layout=True)
contourAx = contourFigure.add_subplot(111, projection=wcs_out)
contourImage = plt.imshow(imageGalactic, cmap='Greys_r', vmax=10)
'''
structures were plotted as contours
each dendrogram level is associated with a color
if dendrogram has more than 6 levels, more colors must be specified
this same affect can be achieved with many levels if you use a color map
this was done previously in commit 564d8e2
'''
colors = ["red", "orange", "yellow", "green",
          "cyan", "blue", "purple", "magenta", "black"]
for i, structure in enumerate(cloudDendrogram):
    if i not in unwantedStructs:
        currentLevel = structure.level
        cloudPlots.plot_contour(contourAx, structure=structure,
                                color=colors[currentLevel])
lon = contourAx.coords[0]
lon.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
lat = contourAx.coords[1]
lat.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
plt.xlim(115, 355)
plt.ylim(115, 355)

# plots beam
my_beam = Beam.from_fits_header(header)
ycen_pix, xcen_pix = 125, 340  # location of beam on plot
pixscale = header['CDELT1']*3600 * u.deg
ellipse_artist = my_beam.ellipse_to_plot(xcen_pix, ycen_pix, pixscale)
plt.gca().add_patch(ellipse_artist)  # plots beam
ellipse_artist.set_facecolor("white")
ellipse_artist.set_edgecolor("black")


# plots scale bar
x = [120, 210]
y = [120, 120]
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

# formats plot
lon.set_ticks(size=-3)
lat.set_ticks(size=-3)
plt.xlabel('Galactic Longtitude', fontsize=10, labelpad=1)
plt.ylabel('Galactic Latitutude', fontsize=10, labelpad=1)
contourAx.tick_params(axis='both', which='major', labelsize=10)
plt.annotate('Continuum', fontsize=10, xy=(0.02, 0.91), xycoords="axes fraction")
plt.savefig("/home/pw/research/Cloud-C/results/continuum/CloudC-structure-contours.pdf")
plt.savefig("/home/pw/research/Cloud-C/results/continuum/CloudC-structure-contours.png")


# continuum and aperture plot
image = fits.getdata(path)
image_2D = np.squeeze(image)
wcs_out, shape_out = find_optimal_celestial_wcs([(image_2D, w1)],
                                                frame='galactic')
imageGalactic, c_footprint = reproject.reproject_interp((image_2D, w1),
                                                        wcs_out,
                                                        shape_out=shape_out)
apertureFigure = plt.figure(2, figsize=(5, 5), constrained_layout=True)
apertureAx = apertureFigure.add_subplot(111, projection=wcs_out)
apertureImage = plt.imshow(imageGalactic, cmap='Greys_r', vmax=5)

lon = apertureAx.coords[0]
lon.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
lat = apertureAx.coords[1]
lat.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
plt.xlim(80, 400)
plt.ylim(80, 400)

# plots beam
my_beam = Beam.from_fits_header(header)
ycen_pix, xcen_pix = 100, 375
pixscale = 0.28 * u.arcsec
ellipse_artist = my_beam.ellipse_to_plot(xcen_pix, ycen_pix, pixscale)
plt.gca().add_patch(ellipse_artist)  # plots beam
ellipse_artist.set_facecolor("white")
ellipse_artist.set_edgecolor("black")

# plots colorbar
cb = plt.colorbar(fraction=0.046, pad=0.04)
cb.set_label(label='Flux Density (mJy / beam)',
             fontsize=10, rotation=270, labelpad=15)
cb.ax.tick_params(which='major', labelsize=10)

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
plt.text(125, 100, '1pc', fontsize=10, color='black')
# formats plot
lon.set_ticks(size=-3)
lat.set_ticks(size=-3)
plt.xlabel('Galactic Longtitude', fontsize=10, labelpad=1)
plt.ylabel('Galactic Latitutude', fontsize=10, labelpad=1)
apertureAx.tick_params(axis='both', which='major', labelsize=10, labelpad=1)
plt.annotate('Continuum', fontsize=10, xy=(0.02, 0.91), xycoords="axes fraction")

# aperture statistics
'''
regions files are just fancy strings, so we create an empty string
that we can append to (remembering to add and endline at the end


I had issues wrt position angle between fk5 and galactic, so I made one regions
file for each frame
'''
regionNames = [1, 7, 2, 4, 5, 3, 6]
print(regionNames)
k = 00
fk5Regions = ''
galRegions = ''
for i, leaf in enumerate(cloudDendrogram.leaves):
    if i not in unwantedLeaves:
        stats = PPStatistic(leaf)
        galacticCoord = wcs_out.pixel_to_world(stats.x_cen, stats.y_cen)
        fk5Coord = galacticCoord.fk5
        # 60.2 degree angular offset between frames (angle of Milky way on sky)
        fk5Region_sky = EllipseSkyRegion(center=fk5Coord,
                                         height=2.3548*stats.minor_sigma.value*pixscale,
                                         width=2.3548*stats.major_sigma.value*pixscale,
                                         angle=stats.position_angle.value*u.deg-60.2*u.deg)
        galacticRegion_sky = EllipseSkyRegion(center=galacticCoord,
                                              height=2.3548*stats.minor_sigma.value*pixscale,
                                              width=2.3548*stats.major_sigma.value*pixscale,
                                              angle=stats.position_angle.value*u.deg)
        newFk5Reg = f'{fk5Region_sky.serialize(format="ds9")}'
        fk5Regions += newFk5Reg[:-2] + f'# text={{Region {regionNames[k]}}} \n'
        newGalReg = f'{galacticRegion_sky.serialize(format="ds9")}'
        galRegions += newGalReg[:-2] + f'# text={{Region {regionNames[k]}}} \n'
        # different from regions files, plots ellipse onto figure
        ellipse = stats.to_mpl_ellipse(edgecolor='red', facecolor='none')
        apertureAx.add_patch(ellipse)
        if k == 2 or k == 3:
            plt.annotate(f"{regionNames[k]}",
                         (stats.x_cen.value, stats.y_cen.value+5))
        if k == 1 or k == 4 or k == 6:
            plt.annotate(f"{regionNames[k]}",
                         (stats.x_cen.value-10, stats.y_cen.value))
        if k == 0 or k == 5:
            plt.annotate(f"{regionNames[k]}",
                         (stats.x_cen.value+5, stats.y_cen.value))
        k += 1
for leaf in cloudCatalogAdjusted:
    region = leaf["_idx"]
    area = leaf["area_ellipse"]
    major = leaf["major_sigma"]
    minor = leaf["minor_sigma"]
    positionAngle = leaf["position_angle"]
    xCen = leaf["x_cen"]
    yCen = leaf["y_cen"]
    print(f"{region}&{area:.3f}&{major:.3f}&{minor:.3f}&{positionAngle:.3f}&{xCen:.3f}&{yCen:.3f}\\\\")
for leaf in cloudCatalogAdjusted:
    region = leaf["_idx"]
    flux = 1e3*leaf["flux"]
    peak = leaf["peak_brightness"]
    print(f"{region}&{flux:.3f}&{peak:.3f}\\\\")
with open("/home/pw/research/Cloud-C/fk5Regions.reg", "w") as table:
    table.write(fk5Regions)
with open("/home/pw/research/Cloud-C/galRegions.reg", "w") as table:
    table.write(galRegions)
plt.savefig("/home/pw/research/Cloud-C/results/continuum/CloudC-structure-apertures.pdf"
            , bbox_inches='tight')
plt.savefig("/home/pw/research/Cloud-C/results/continuum/CloudC-structure-apertures.png"
            , bbox_inches='tight')
