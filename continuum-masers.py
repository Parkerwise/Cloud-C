'''
-*- coding: utf-8 -*-
Author: Parker Wise
Date: 2024-06-27
Description: Plots masers onto continuum
Python Version: 3.11.9
'''

import astropy.io.fits as fits  # 6.1.0
from astropy.coordinates import SkyCoord
import pandas as pd  # 2.2.2
import matplotlib.pyplot as plt  # 3.8.0
import matplotlib.colors as colors
from astropy.wcs import WCS
from radio_beam import Beam  # 0.3.7
from astropy import units as u
import numpy as np  # 1.26.4
import reproject  # 0.13.1
from reproject.mosaicking import find_optimal_celestial_wcs
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
tick_font_size = 10
plt.rcParams['text.usetex'] = True
path = "/home/pw/research/Cloud-C/co-data/CloudC_mJy.fits"
image = fits.getdata(path)
header = fits.getheader(path)
w1 = WCS(header)
w1 = w1.dropaxis(3)
w1 = w1.dropaxis(2)
image_2D = np.squeeze(image)
# plots in galactic
wcs_out, shape_out = find_optimal_celestial_wcs([(image_2D, w1)],
                                                frame='galactic')
cont, c_footprint = reproject.reproject_interp((image_2D, w1), wcs_out,
                                               shape_out=shape_out)

# init plot and axes
fig1 = plt.figure(1, figsize=(10, 10), constrained_layout=True)
ax1 = plt.subplot(projection=wcs_out)
lon = ax1.coords[0]
lon.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
lat = ax1.coords[1]
lat.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
plt.xlim(80, 410)
plt.ylim(80, 410)

# formats beam
im1 = plt.imshow(cont, cmap='Greys_r', vmax=5)  # plots continuum
my_beam = Beam.from_fits_header(header)
ycen_pix, xcen_pix = 100, 375
pixscale = 0.28 * u.arcsec
ellipse_artist = my_beam.ellipse_to_plot(xcen_pix, ycen_pix, pixscale)
plt.gca().add_patch(ellipse_artist)  # plots beam
ellipse_artist.set_facecolor("white")
ellipse_artist.set_edgecolor("black")


# plots scalebar
x = [90, 180]
y = [90, 90]
plt.plot(np.array([x[0], x[1]]), np.array([y[0], y[1]]),
         color="black", linewidth=3)
scalebarBegin = w1.pixel_to_world(x[0], y[0])
scalebarEnd = w1.pixel_to_world(x[1], y[1])
sep = scalebarBegin.separation(scalebarEnd)
plt.text(125, 100, '1pc', fontsize=14, color='black')

# formats plot
lon.set_ticks(size=-3)
lat.set_ticks(size=-3)
plt.xlabel('Galactic Longtitude', fontsize=20, labelpad=1)
plt.ylabel('Galactic Latitutude', fontsize=20, labelpad=1)
ax1.tick_params(axis='both', which='major', labelsize=15)
plt.annotate('Continuum', fontsize=15, xy=(0.02, 0.91),
             xycoords="axes fraction")

# masers
# header skips over comments
df = pd.read_csv('/home/pw/research/Cloud-C/co-data/masers.csv', header=2)
pixScale = 0.28  # asec per pixel


def asec2pix(asec):
    return asec/pixScale


def coord2pixel(lon, lat):
    skycoord = SkyCoord(lon, lat, unit="deg", frame="galactic")
    pixelcoord = skycoord.to_pixel(wcs_out, 0, mode="wcs")
    return pixelcoord[0], pixelcoord[1]  # returns x, y positions


positions = [coord2pixel(df.l[i], df.b[i]) for i in range(len(df.l))]
l_err = [asec2pix(sigma_l) for sigma_l in df.sigma_l]
b_err = [asec2pix(sigma_b) for sigma_b in df.sigma_b]


x = [positions[i][0] for i in range(len(positions))]
y = [positions[i][1] for i in range(len(positions))]

cmap = plt.cm.jet
# markers were made to reflect the markers used in Ginsburg 2015
markers = ["p", "o", "^", "v", "D", "s", "s", "s", "s", "s", "s", "s", "*"]
names = ['CH$_3$OH $7_0 - 6_1\\mathrm{ A}^+$',
         'H$_2$CO $1_{1,0}-1_{1,1}$',
         'SiO J=1-0 v=1',
         'SiO J=1-0 v=2',
         'CH$_3$OH $5_1 - 6_0\\mathrm{ A}^+$',
         'H$_2$O',
         '',
         '',
         '',
         '',
         '',
         '',
         'OH'
         ]
# to use scatter with intensity we have to normalize our data
# vmin --> 0, vmax --> 1
norm = colors.Normalize(vmin=9, vmax=80)
for x, y, vel, mark, l_err, b_err, name in zip(x, y, df.velocity,
                                               markers, l_err, b_err, names):
    scatter = plt.scatter(x=x, y=y, c=cmap(norm(vel)), s=80,
                          marker=mark, zorder=1, alpha=0.75, label=name)
    plt.errorbar(x=x, y=y, yerr=b_err, xerr=l_err, fmt="o", zorder=0,
                 color="black", lw=3)
# ScalarMappable is needed to scale the color bar correctly
scatter = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
scatter._A = np.array([norm.vmin, norm.vmax])
scatterBar = plt.colorbar(scatter, fraction=0.046, pad=0.04, ax=plt.gca())
scatterBar.set_label('Velocity (km/s)', fontsize=25, rotation=270, labelpad=30)
plt.legend()
# always save pdf and png! pdf work well in papers
# but sometimes a png is the better option (in slideshow)
plt.savefig("/home/pw/research/Cloud-C/results/continuum/continuum-masers.pdf",
            dpi=250, pad_inches=1)
plt.savefig("/home/pw/research/Cloud-C/results/continuum/continuum-masers.png",
            dpi=250, pad_inches=1)
plt.show()
