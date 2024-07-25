'''
-*- coding: utf-8 -*-
Author: Parker Wise
Date: 2024-07-05
Description: Masks continuum sources and calculates standard deviation
Python Version: 3.11.9
'''

import astropy.io.fits as fits  # 6.1.0
from astropy.wcs import WCS
import numpy as np  # 1.26.4
import pandas as pd  # 2.2.2
import ast  # 0.9.33
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
# standard way of importing continuum as ndarray
path = "/home/pw/research/Cloud-C/co-data/CloudC_mJy.fits"
image = fits.getdata(path)
header = fits.getheader(path)
w1 = WCS(header)
w1 = w1.dropaxis(3)
w1 = w1.dropaxis(2)
image_2D = np.squeeze(image)

'''
this code does not need to import a CSV, if you manually call
addMask it'll still work. I'm just importing a CSV since we've ID'd the
sources for this cloud by eye
'''
sourcesPath = "/home/pw/research/CO-ACES/Tables/cloud_c_locations.csv"
df = pd.read_csv(sourcesPath, usecols=[1, 4], names=["center", "radius"],
                 header=1)
# pandas imports strings. ast takes strings --> tuples
centers = [ast.literal_eval(entry) for entry in df.center]


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


for center, radius in zip(centers, df.radius):
    addMask(image_2D, radius, center)

addMask(image_2D, 100, (166, 175), remove_inside=False)  # crops noisy edges
sigma = np.nanstd(image_2D)
'''
# for viewing mask
import matplotlib.pyplot as plt
im1 = plt.imshow(image_2D, cmap='Greys_r', vmax=5)  # plots continuum
plt.show()
'''
print(f"sigma: {sigma}, 5 sigma: {5*sigma}")
