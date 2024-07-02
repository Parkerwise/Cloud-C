# -*- coding: utf-8 -*- 
# Author: Parker Wise
# Date: 2024-07-01
# Description: Creates CSV of spectral cube properties from several files
# Python Version: 3.11.9
import astropy.io.fits as fits #6.1.0
from decimal import Decimal
import numpy as np #1.26.4

fitsFilenames=["A.Dust_Ridge_12C18O.cube.I.pbcor.fits",
               "A.Dust_Ridge_13C16O.cube.I.pbcor.fits",
               "A.Dust_Ridge_sci.spw29.cube.I.pbcor.fits",
               "A.Dust_Ridge_sci.spw31.cube.I.pbcor.fits",
               "B.Dust_Ridge_12C16O_1-0.cube.I.pbcor.fits",
               "B.Dust_Ridge_12C17O.cube.I.pbcor.fits",
               "B.Dust_Ridge_sci.spw25_27_29_31.cont.I.tt0.pbcor.fits",
               "B.Dust_Ridge_sci.spw29.cube.I.pbcor.fits",
               "B.Dust_Ridge_sci.spw31.cube.I.pbcor.fits"]
#allows us to find header key given header value
def getHeaderList(header):
    try:
        headerKeys=np.array(list(header.keys()))
        headerValues=np.array(list(header.values()))
    except AttributeError:
        print("AttributeError: FITS file must have a header")
    return headerKeys, headerValues

def fetchFreqAxis(header): #finds what axis freq is
    headerKeys, headerValues=getHeaderList(header)
    #it is assumed there is only one frequency axis
    freqIndex=np.where(#or operators allows for all spectral types
            (headerValues=='FREQ') | 
            (headerValues=='ENER') |
            (headerValues=='WAVN') |
            (headerValues=='VRAD') |
            (headerValues=='WAVE') |
            (headerValues=='VOPT') |
            (headerValues=='ZOPT') |
            (headerValues=='AWAV') |
            (headerValues=='VELO') |
            (headerValues=='BETA') 
    )
    #returns false if there is no spectral axis
    freqKey=headerKeys[freqIndex]
    if len(freqKey)==0: 
        return False 
    else:
        freqAxes=freqKey[0][-1]
    return freqAxes

def fetchSpatialAxes(header):
    headerKeys, headerValues=getHeaderList(header)
    spatialIndices=np.where(headerValues=='deg')
    spatialKeys=headerKeys[spatialIndices]
    if len(spatialKeys)==0:
        return False 
    else:
        spatialAxes=[key[-1] for key in spatialKeys]
    return spatialAxes

def cube_properties(path):
    image=fits.getdata(path) 
    header=fits.getheader(path)
    calibration=path[0]
    spectralWindow=header['FILNAM04']
    name=calibration+'.'+spectralWindow
    xdim=header['CTYPE1']+": "+header['CUNIT1']
    ydim=header['CTYPE2']+": "+header['CUNIT2']
    zdim=header['CTYPE3']+": "+header['CUNIT3']
    spatialAxes=fetchSpatialAxes(header)
    if spatialAxes==False: 
        pixelScale="---"
    else:
        #pixel scale is multiplied from degrees to asec
        pixelScale=f"{Decimal(abs(header['CDELT1'])*3600):.3E} asec/px"
    bmaj=f"{Decimal(header['BMAJ']*3600):.3E} asec/px"
    bmin=f"{Decimal(header['BMIN']*3600):.3E} asec/px"
    bpa=f"{Decimal(header['BPA']):.3E} deg"

    freqAxis=fetchFreqAxis(header)
    if freqAxis==False:
        freqRange="---"
        freqWidth="---"
        channels="---"
        channelWidth="---"
    else:
        freqUnit=header['CUNIT'+freqAxis]
        #image is ordered in reverse: CTYPE4,CTYPE3...
        channels=np.shape(image)[-int(freqAxis)] 
        #frequencies are converted to GHz instead of Hz
        startFreqVal=header['CRVAL'+freqAxis]*10**-9
        channelWidth=header['CDELT'+freqAxis]*10**-9
        endFreq=startFreqVal+channelWidth*channels
        freqWidth=(endFreq-startFreqVal)
        #Decimal formats values in scientific notation
        freqWidth=f"{Decimal(freqWidth):.3E} {freqUnit}"
        startFreq=f"{Decimal(float(startFreqVal)):.3E}"
        endFreq=f"{Decimal(float(endFreq)):.3E}"
        freqRange=f"[{startFreq} - {endFreq}] G{freqUnit}"
        channelWidth=f"{Decimal(channelWidth):.3E} G{freqUnit}"
    return startFreqVal, f"{name},{xdim}, {ydim},{zdim}, {pixelScale}, {bmaj}, {bmin}, {bpa}, {channels},{channelWidth},{freqWidth}, {freqRange}\n"
#csv should be sorted by ascending frequency
#to avoid having to open all of our files twice, we save the freq and the
#desired output to a dictionary that we then sort by frequency
entryAndFreq={}
for path in fitsFilenames:
    freq, entry=cube_properties(path)
    entryAndFreq[freq]=entry
sortedEntries=sorted(entryAndFreq.items())
with open("cube-properties.csv",'w') as table: 
    table.write("name, xdim, ydim, zdim, pixelScale, bmaj, bmin, bpa, channels,channelWidth,freqWidth, freqRange\n")
    for entry in sortedEntries:
        table.write(entry[1])
