'''
-*- coding: utf-8 -*-
Author: Parker Wise
Date:
Description:
Python Version: 3.11.9
'''
from astropy.wcs import WCS
import matplotlib.colors as colors
from radio_beam import Beam
from regions import Regions
import reproject
from reproject.mosaicking import find_optimal_celestial_wcs
import astropy.io.fits as fits
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
import sys
import matplotlib
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# Reading in files
# # This should be the file where we're extracting spectra
matplotlib.rcParams.update({'font.size': 10})
data_dir = "/home/pw/research/Cloud-C/co-data/"

cube_list = [
    "A.Dust_Ridge_12C18O.cube.I.pbcor.fits",
    "A.Dust_Ridge_13C16O.cube.I.pbcor.fits",
    "B.Dust_Ridge_12C16O_1-0.cube.I.pbcor.fits",
    "B.Dust_Ridge_12C17O.cube.I.pbcor.fits",
]

cube_abbreviation = [
    "ASPW25",  # 12C18O
    "ASPW27",  # 13C16O
    "BSPW25",  # 12C16O
    "BSPW27",  # 12C17O
]
for path, name in zip(cube_list, cube_abbreviation):
    path = data_dir + path
    print(name)
    header = fits.getheader(path)
    w1 = WCS(header)
    w1 = w1.dropaxis(3)
    w1 = w1.dropaxis(2)
    # # creates spectral cube object
    sc = SpectralCube.read(path)
    sc.allow_huge_operations = True

    # Converts units to GHz and K, feel free to change
    # #  defaults should be Hz and Jy/Beam
    sc_kms = sc.with_spectral_unit(u.km/u.s, velocity_convention="radio")
    print(sc_kms)
    sc_slab = sc_kms.spectral_slab(-100. * u.km / u.s, 200. * u.km / u.s)
    sc_slab.allow_huge_operations = True
    sc_bin = sc_slab.downsample_axis(3, axis=0)
    sc_K_kms = sc_bin.to(u.K)

    moment_0 = sc_K_kms.moment(order=0, how='slice')
    moment_1 = sc_K_kms.moment(order=1, how='slice')
    moment_0_max = np.nanmax(moment_0.hdu.data)
    moment_0_std = np.nanstd(moment_0.hdu.data)
    threshold = 1.25 * moment_0_std
    badpix = np.where(moment_0.hdu.data < threshold)
    moment_1.hdu.data[badpix] = np.nan
    fig1 = plt.figure(1, figsize=(5, 5), constrained_layout=True)

    wcs_out, shape_out = find_optimal_celestial_wcs([(moment_1.hdu.data, w1)],
                                                    frame='galactic')
    # Rotate image
    moment_0_cont, moment_0_c_footprint = reproject.reproject_interp((moment_0.hdu.data, w1),
                                                                      wcs_out, shape_out=shape_out)
    cont, c_footprint = reproject.reproject_interp((moment_1.hdu.data, w1),
                                                   wcs_out, shape_out=shape_out)

    ax1 = plt.subplot(1, 1, 1, projection=wcs_out)
    lon = ax1.coords[0]
    lon.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
    lat = ax1.coords[1]
    lat.set_format_unit(u.deg, decimal=True, show_decimal_unit=True)
    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    padding = 20
    bound = np.argwhere(~np.isnan(moment_0_cont))
    x_min = min(bound[:, 1]-padding)
    y_min = min(bound[:, 0]-padding)
    x_max = max(bound[:, 1]+padding)
    y_max = max(bound[:, 0]+padding)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)

    my_beam = Beam.from_fits_header(header)
    ycen_pix, xcen_pix = y_min+20, x_max-20
    pixscale = 0.28 * u.arcsec
    ellipse_artist = my_beam.ellipse_to_plot(xcen_pix, ycen_pix, pixscale)
    mean_velocity = np.nanmean(cont)
    image_std = np.nanstd(cont)
    # vmin = mean_velocity - 2 * image_std
    # vmax = mean_velocity + 2 * image_std
    vmin = -25
    vmax = 100
    print(vmin, vmax)
    im1 = plt.imshow(cont, cmap='rainbow', vmin=vmin, vmax=vmax)
    plt.gca().add_patch(ellipse_artist)
    ellipse_artist.set_facecolor("white")
    ellipse_artist.set_edgecolor("black")

    cmap = plt.cm.rainbow
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    scatter = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    scatter._A = np.array([norm.vmin, norm.vmax])
    cb = plt.colorbar(scatter, fraction=0.046, pad=0.04, ax=plt.gca())
    cb.set_label(label='Mean Weighted Velocity (km/s)',
                 fontsize=10, rotation=270, labelpad=10)
    cb.ax.tick_params(which='major', labelsize=10, pad=10)


    plt.annotate(f'Moment 1\n{name}', fontsize=10, xy=(0.02, 0.91),
                 xycoords="axes fraction")
    x = [x_min+10, x_min+100]
    y = [y_min+10, y_min+10]
    plt.plot(np.array([x[0], x[1]]),
             np.array([y[0], y[1]]),
             color="black", linewidth=3)
    scalebarBegin = w1.pixel_to_world(x[0], y[0])
    scalebarEnd = w1.pixel_to_world(x[1], y[1])
    plt.xlabel('Galactic Longtitude', fontsize=10, labelpad=1)
    plt.ylabel('Galactic Latitutude', fontsize=10, labelpad=1)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    plt.text(x_min+45, y_min+20, '1pc', fontsize=10, color='black')

    regions_file = "/home/pw/research/Cloud-C/galRegions.reg"
    regionNames = [1, 7, 2, 4, 5, 3, 6]
    regs = Regions.read(regions_file)
    k = 0
    for sky_region in regs:
        pixel_region = sky_region.to_pixel(wcs_out)
        pixel_region.plot(color="black")
        x = pixel_region.center.xy[0]
        y = pixel_region.center.xy[1]
        if k == 2 or k == 3:
            plt.annotate(f"{regionNames[k]}",
                         (x-5, y-20))
        if k == 1 or k == 4 or k == 6:
            plt.annotate(f"{regionNames[k]}",
                         (x-15, y))
        if k == 0 or k == 5:
            plt.annotate(f"{regionNames[k]}",
                         (x+5, y))
        k += 1

    plt.savefig(f"results/moment-maps/{name}-moment_1.pdf")
    plt.savefig(f"results/moment-maps/{name}-moment_1.png", dpi=200)
