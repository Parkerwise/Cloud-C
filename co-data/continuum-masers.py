import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
import pandas as pd
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.pyplot as plt    
import matplotlib.colors as colors
from astropy.wcs import WCS    
from radio_beam import Beam
from astropy import units as u  
import numpy as np
import reproject
from reproject.mosaicking import find_optimal_celestial_wcs 
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

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

#plots beam
my_beam = Beam.from_fits_header(header)  
ycen_pix, xcen_pix = 100, 375
pixscale = 0.28 * u.arcsec
ellipse_artist = my_beam.ellipse_to_plot(xcen_pix, ycen_pix, pixscale)
im1 = plt.imshow(cont,cmap='Greys_r',vmax=5)
plt.gca().add_patch(ellipse_artist)
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

#masers
df = pd.read_csv('masers.csv',header=2) #header skips over comments
dist=8.2*10**3 
pixScale=0.28 # asec per pixel
def asec2pix(asec):
    return asec/pixScale
def coord2pixel(l,b):
    skycoord = SkyCoord(l, b, unit="deg", frame="galactic")
    pixelcoord=skycoord.to_pixel(wcs_out,0,mode="wcs")
    x_pixel=pixelcoord[0]
    y_pixel=pixelcoord[1]
    return x_pixel,y_pixel

positions=[coord2pixel(df.l[i],df.b[i]) for i in range(len(df.l))]
l_err=[asec2pix(sigma_l) for sigma_l in df.sigma_l]
b_err=[asec2pix(sigma_b) for sigma_b in df.sigma_b]
vel=np.asarray([39.6, 36.7, 64.8, 63.2, 38.0, 78.8, 9.5, 24.4, 32.3, 37.3, 40.4, 52.5, 36.0])
x=[positions[i][0] for i in range(len(positions))]
y=[positions[i][1] for i in range(len(positions))]
cmap=plt.cm.jet
markers=["p", "o", "^", "v", "D", "s","s","s","s","s","s","s", "*"]
norm = colors.Normalize(vmin=9, vmax=80)
for x,y,velocity,mark,l_err,b_err in zip(x,y,vel,markers,l_err,b_err):
    print(velocity)
    sc = plt.scatter(x=x,y=y,c=cmap(norm(velocity)),s = 80 ,marker = mark,zorder=1,alpha=0.75)
    plt.errorbar(x=x,y=y,yerr=b_err, xerr=l_err, fmt="o",zorder=0,color="black",lw=3)
sc = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sc._A = np.array([norm.vmin, norm.vmax])
cb2=plt.colorbar(sc,fraction=0.046,pad=0.04,ax=plt.gca())                                      
cb2.set_label('Velocity (km/s)',fontsize=25,rotation=270,labelpad=30)
cb=plt.colorbar(im1,fraction=0.046,pad=0.04)                                      
cb.set_label(label='Flux Density (mJy / beam)',fontsize=25,rotation=270,labelpad=30) 
cb.ax.tick_params(which = 'major', labelsize = 20) 
plt.legend()
plt.show()
#saves fig
#plt.savefig("continuum.pdf",dpi=250,pad_inches=1)
#plt.savefig("continuum.png",dpi=250,pad_inches=1)
