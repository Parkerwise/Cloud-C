#CODE IS A WORK IN PROGRESS

# -*- coding: utf-8 -*- 
# Author: Parker Wise
# Date: 
# Description: 
# Python Version: 3.11.9

#librarys
from astrodendro import Dendrogram, pp_catalog
import astropy.io.fits as fits #6.1.0
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt #3.6.0
import matplotlib.colors as colors
from astropy.wcs import WCS    
from radio_beam import Beam #0.3.7
from astropy import units as u  
import numpy as np #1.26.4
import reproject #0.13.1
from reproject.mosaicking import find_optimal_celestial_wcs 
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
#plotting and stuff
#formatting stuff
plt.rcParams['text.usetex'] = True
path ="CloudC_mJy.fits"
image=fits.getdata(path)
header=fits.getheader(path)
w1=WCS(header)
w1 = w1.dropaxis(3)
w1 = w1.dropaxis(2)
image_2D = np.squeeze(image)
#plots in galactic
wcs_out,shape_out= find_optimal_celestial_wcs([(image_2D,w1)],frame='galactic')
cont,c_footprint = reproject.reproject_interp((image_2D,w1),wcs_out,shape_out=shape_out)
 
#init plot and axes
'''
fig1 = plt.figure(1,figsize=(10,10),constrained_layout=True)       
ax1 = plt.subplot(projection=wcs_out)             
lon = ax1.coords[0]
lon.set_format_unit(u.deg,decimal=True,show_decimal_unit=True)
lat = ax1.coords[1]
lat.set_format_unit(u.deg,decimal=True,show_decimal_unit=True)
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
plt.xlim(80,410)
plt.ylim(80,410)

#formats beam
my_beam = Beam.from_fits_header(header)  
ycen_pix, xcen_pix = 100, 375
pixscale = 0.28 * u.arcsec
ellipse_artist = my_beam.ellipse_to_plot(xcen_pix, ycen_pix, pixscale)
im1 = plt.imshow(cont,cmap='Greys_r',vmax=5) #plots continuum
plt.gca().add_patch(ellipse_artist) #plots beam
ellipse_artist.set_facecolor("white")
ellipse_artist.set_edgecolor("black")


#plots scalebar
x=[90,180]
y=[90,90]
plt.plot(np.array([x[0],x[1]]),np.array([y[0],y[1]]),color="black",linewidth=3)
scalebarBegin=w1.pixel_to_world(x[0],y[0])
scalebarEnd=w1.pixel_to_world(x[1],y[1])
sep=scalebarBegin.separation(scalebarEnd)
plt.text(125,100,'1pc',fontsize=14,color='black')

#formats plot
lon.set_ticks(size=-3)                                                                                     
lat.set_ticks(size=-3)                                                                                   
plt.xlabel('Galactic Longtitude',fontsize=20,labelpad=1)                              
plt.ylabel('Galactic Latitutude',fontsize=20,labelpad=1)                                  
ax1.tick_params(axis = 'both', which = 'major', labelsize = 15)                      
plt.annotate('Continuum',fontsize=15,xy=(0.02,0.91),xycoords="axes fraction")
'''

def addMask(image,radius,center,remove_inside=True):
    Y,X=np.ogrid[:350,:350] 
    dist=np.sqrt((X-center[0])**2+(Y-center[1])**2) 
    mask= dist <= radius
    if remove_inside==True:
        for (i,j),bool in np.ndenumerate(mask):
            if bool == True:
                image[i,j]=np.nan
    if remove_inside==False:
        for (i,j),bool in np.ndenumerate(mask):
            if bool == False:
                image[i,j]=np.nan

addMask(image_2D,110,(166,175),remove_inside=False) #crops noisy edges
sigma=0.41 #mJy/beam
#pixel per beam (do later) 
d=Dendrogram.compute(image_2D,min_value=sigma,min_delta=2*sigma)
metadata={ }

#plt.savefig("continuum-masers.pdf",dpi=250,pad_inches=1)
#plt.savefig("continuum-masers.png",dpi=250,pad_inches=1)



