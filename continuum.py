import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from radio_beam import Beam
import astropy.units as u
import pylab
import numpy as np
import reproject
from reproject.mosaicking import find_optimal_celestial_wcs
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
'''
path ="B.Dust_Ridge_sci.spw25_27_29_31.cont.I.tt0.pbcor.fits"
image=fits.getdata(path)
mJy_image = image * 1000        # 1 Jy = 1000 mJy so we would multiply our image by 1000
header['BUNIT'] = 'mJy/beam '
print(header['BUNIT'])
fits.writeto('CloudC_mJy.fits',mJy_image,header,overwrite=True)
'''

plt.rcParams['xtick.major.pad'] = '8'
path = "co-data/CloudC_mJy.fits"
image = fits.getdata(path)

header = fits.getheader(path)
w1 = WCS(header)
w1 = w1.dropaxis(3)
w1 = w1.dropaxis(2)

image_2D = pylab.squeeze(image)

######
# Use reproject functions
wcs_out, shape_out = find_optimal_celestial_wcs([(image_2D, w1)],
                                                frame='galactic')
# Rotate image
cont, c_footprint = reproject.reproject_interp((image_2D, w1),
                                               wcs_out, shape_out=shape_out)
# init plot and axes
fig1 = plt.figure(1, figsize=(5, 5), constrained_layout=True)
ax1 = plt.subplot(projection=wcs_out)
lon = ax1.coords[0]
lon.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
lat = ax1.coords[1]
lat.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
plt.xlim(80, 400)
plt.ylim(80, 400)
my_beam = Beam.from_fits_header(header)
ycen_pix, xcen_pix = 100, 375
pixscale = 0.28 * u.arcsec
ellipse_artist = my_beam.ellipse_to_plot(xcen_pix, ycen_pix, pixscale)
im1 = plt.imshow(cont, cmap='Greys_r', vmax=5)
plt.gca().add_patch(ellipse_artist)
ellipse_artist.set_facecolor("white")
ellipse_artist.set_edgecolor("black")

x = [90, 180]
y = [90, 90]
plt.plot(np.array([x[0], x[1]]),
         np.array([y[0], y[1]]),
         color="black", linewidth=3)
scalebarBegin = w1.pixel_to_world(x[0], y[0])
scalebarEnd = w1.pixel_to_world(x[1], y[1])
sep = scalebarBegin.separation(scalebarEnd)
print(sep.degree)
print(sep.degree*3600)
pylab.text(125, 100, '1pc', fontsize=10, color='black')

# formats plot
pylab.xlabel('Galactic Longtitude', fontsize=10, labelpad=1)
pylab.ylabel('Galactic Latitutude', fontsize=10, labelpad=1)
ax1.tick_params(axis='both', which='major', labelsize=10)
plt.annotate('Continuum', fontsize=10, xy=(0.02, 0.91),
             xycoords="axes fraction")
cb = pylab.colorbar(im1, fraction=0.046, pad=0.04)
cb.set_label(label='Flux Density (mJy / beam)',
             fontsize=10, rotation=270, labelpad=10)
cb.ax.tick_params(which='major', labelsize=10, pad=10)

# saves fig
pylab.savefig("results/continuum/continuum.pdf", dpi=250, pad_inches=1)
pylab.savefig("results/continuum/continuum.png", dpi=250, pad_inches=1)
