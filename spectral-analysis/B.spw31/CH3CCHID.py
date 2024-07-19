#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# Author: Parker Wise
# Date: 2024-06-21
# Description: fits a model to several transitions of CH3CCH in cloud 'c'
# Python Version: 3.12.4
import astropy.io.fits as fits  #7.0.0
import matplotlib.pyplot as plt #3.9.0   
from astropy.wcs import WCS
from astropy import units as u  
import pyspeckit as psk #1.0.4
from regions import Regions #0.9
import numpy as np #2.0.0
from spectral_cube import SpectralCube #0.6.5
import sys
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

path="/home/pw/research/Cloud-C/co-data/B.Dust_Ridge_sci.spw31.cube.I.pbcor.fits"
#spectrum is calculated for the area inside region file
regions_file="/home/pw/research/Cloud-C/co-data/c1-by-eye" 
header=fits.getheader(path)
w1=WCS(header)

#obtaining frequencies and brightnesses
sc=SpectralCube.read(path)
sc.allow_huge_operations=True 
sc_Ghz=sc.with_spectral_unit(u.GHz)
sc_Ghz=sc_Ghz.to(u.K)
freq,Dec,Ra = sc_Ghz.world[:,0,0]
regpix = Regions.read(regions_file)
subcube = sc_Ghz.subcube_from_regions([regpix[0]])  
spectrum = subcube.mean(axis=(1, 2))

#the STD is needed for the fitting
Noise_upper=825 #defines some line free channels to calculate STD
Noise_lower=800
sigma=np.std(spectrum[Noise_lower:Noise_upper].value)

freq_min,freq_max=102.45,102.57 #the chunk of spectrum we're analyzing
zoom = np.where((freq.value>freq_min)*(freq.value<freq_max))
zoom_spectrum=spectrum[zoom]
zoom_freq = freq[zoom]
error = np.zeros(np.size(zoom))
meas = sigma #uses sigma calc from previous block
error = error + meas
sp=psk.Spectrum(data=zoom_spectrum,xarr=zoom_freq,error=error,unit='K')

#line rich regions to exclude when calculating baseline
exclude=np.array([102.47,102.56])

#heights and line widths were guessed by eye
#freqs were calculated for the transitions of CH3CCH at 40 km/s
         #Height  #center freq    #line width
guesses=[2,       102.534,        0.001, #6(0)-5(0)
         1,       102.532,        0.001, #6(1)-5(1)
         1,       102.527,        0.001, #6(2)-5(2)
         1,       102.517,        0.001, #6(3)-5(3)
         0.25,    102.503,        0.001] #6(4)-5(4)

#limits give acceptable ranges of variation
        #Height     #center freq       #line width
limits=[(1.50,3.0), (102.525,102.540), (0.0001,0.01),  #6(0)-5(0)
        (0.75,1.5), (102.525,102.535), (0.0001,0.01),  #6(1)-5(1)
        (0.75,1.5), (102.520,102.530), (0.0001,0.01),  #6(2)-5(2)
        (0.75,1.5), (102.510,102.525), (0.0001,0.01),  #6(3)-5(3)
        (0.25,1.0), (102.500,102.510), (0.0001,0.01)]  #6(4)-5(4)

#false limited value means to ignore limits
limited = [(True,True),(True,True),(True,True)
           ,(True,True),(True,True),(True,True)
           ,(True,True),(True,True),(True,True)
           ,(True,True),(True,True),(True,True)
           ,(True,True),(True,True),(True,True)]
sp.baseline(order=1,exclude=exclude)#Fit and subtract a polynomial (order=1 -> linear) baseline
#sp.specfit(guesses=guesses,limits=limits,limited=limited)
sp.specfit(guesses=guesses)


#fig 1 shows the fit of the model to the data
fig1 = plt.figure(1,figsize=(10,5),dpi=250)
sp.plotter(axis=plt.subplot(1,1,1),errstyle='fill')
sp.specfit.plot_fit(axis=plt.subplot(1,1,1),annotate=False)
sp.specfit.plotresiduals(axis=plt.subplot(1,1,1),clear=False, yoffset=-2,label=False)
#to identify the transitions within figures
plt.vlines(sp.specfit.parinfo[1].value,0,sp.specfit.parinfo[0].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[4].value,0,sp.specfit.parinfo[3].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[7].value,0,sp.specfit.parinfo[6].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[10].value,0,sp.specfit.parinfo[9].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[13].value,0,sp.specfit.parinfo[12].value+0.3,colors="green",label="",linewidth=1)
plt.text(sp.specfit.parinfo[1].value-0.001,sp.specfit.parinfo[0].value+0.4,'6(0)-5(0)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[4].value-0.001,sp.specfit.parinfo[3].value+0.4,'6(1)-5(1)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[7].value-0.001,sp.specfit.parinfo[6].value+0.4,'6(2)-5(2)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[10].value-0.001,sp.specfit.parinfo[9].value+0.4,'6(3)-5(3)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[13].value-0.001,sp.specfit.parinfo[12].value+0.4,'6(4)-5(4)',rotation='vertical',color='green')
plt.ylim(-0.1,3)
sp.plotter.savefig('CH3CCH.fit.png')
sp.plotter.savefig('CH3CCH.fit.pdf')

#fig 2 shows each component
fig2 = plt.figure(2,figsize=(10,5),dpi=250)
sp.plotter(axis=plt.subplot(1,1,1),errstyle='fill')
sp.specfit.plot_components(axis=plt.subplot(1,1,1),add_baseline=False,component_yoffset=0)
sp.specfit.plotresiduals(axis=plt.subplot(1,1,1),clear=False, yoffset=-2,label=False)
plt.vlines(sp.specfit.parinfo[1].value,0,sp.specfit.parinfo[0].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[4].value,0,sp.specfit.parinfo[3].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[7].value,0,sp.specfit.parinfo[6].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[10].value,0,sp.specfit.parinfo[9].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[13].value,0,sp.specfit.parinfo[12].value+0.3,colors="green",label="",linewidth=1)
plt.text(sp.specfit.parinfo[1].value-0.001,sp.specfit.parinfo[0].value+0.4,'6(0)-5(0)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[4].value-0.001,sp.specfit.parinfo[3].value+0.4,'6(1)-5(1)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[7].value-0.001,sp.specfit.parinfo[6].value+0.4,'6(2)-5(2)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[10].value-0.001,sp.specfit.parinfo[9].value+0.4,'6(3)-5(3)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[13].value-0.001,sp.specfit.parinfo[12].value+0.4,'6(4)-5(4)',rotation='vertical',color='green')
plt.ylim(-0.1,3)
sp.plotter.savefig('CH3CCH.components.png')
sp.plotter.savefig('CH3CCH.components.pdf')

#fig 3 shows both the fit and the components
fig3 = plt.figure(3,figsize=(10,5),dpi=250)
sp.plotter(axis=plt.subplot(1,1,1),errstyle='fill')
sp.specfit.plot_fit(axis=plt.subplot(1,1,1),annotate=False)
sp.specfit.plot_components(axis=plt.subplot(1,1,1),add_baseline=False,component_yoffset=-0.3)
sp.specfit.plotresiduals(axis=plt.subplot(1,1,1),clear=False, yoffset=-2,label=False)
plt.vlines(sp.specfit.parinfo[1].value,0,sp.specfit.parinfo[0].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[4].value,0,sp.specfit.parinfo[3].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[7].value,0,sp.specfit.parinfo[6].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[10].value,0,sp.specfit.parinfo[9].value+0.3,colors="green",label="",linewidth=1)
plt.vlines(sp.specfit.parinfo[13].value,0,sp.specfit.parinfo[12].value+0.3,colors="green",label="",linewidth=1)
plt.text(sp.specfit.parinfo[1].value-0.001,sp.specfit.parinfo[0].value+0.4,'6(0)-5(0)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[4].value-0.001,sp.specfit.parinfo[3].value+0.4,'6(1)-5(1)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[7].value-0.001,sp.specfit.parinfo[6].value+0.4,'6(2)-5(2)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[10].value-0.001,sp.specfit.parinfo[9].value+0.4,'6(3)-5(3)',rotation='vertical',color='green')
plt.text(sp.specfit.parinfo[13].value-0.001,sp.specfit.parinfo[12].value+0.4,'6(4)-5(4)',rotation='vertical',color='green')
plt.ylim(-0.5,3)
sp.plotter.savefig('CH3CCH.png')
sp.plotter.savefig('CH3CCH.pdf')

#so values can be looked at later on
with open("CH3CCH.transitions.csv", "w") as table:
    table.write('Transition,Central Freq (GHz),Amplitude (K),Width (GHz)\n')
    table.write(f'6(0)-5(0),{sp.specfit.parinfo[1].value},{sp.specfit.parinfo[0].value},{sp.specfit.parinfo[2].value}\n')
    table.write(f'6(1)-5(1),{sp.specfit.parinfo[4].value},{sp.specfit.parinfo[3].value},{sp.specfit.parinfo[5].value}\n')
    table.write(f'6(2)-5(2),{sp.specfit.parinfo[7].value},{sp.specfit.parinfo[6].value},{sp.specfit.parinfo[8].value}\n')
    table.write(f'6(3)-5(3),{sp.specfit.parinfo[10].value},{sp.specfit.parinfo[9].value},{sp.specfit.parinfo[11].value}\n')
    table.write(f'6(4)-5(4),{sp.specfit.parinfo[13].value},{sp.specfit.parinfo[12].value},{sp.specfit.parinfo[14].value}\n')
