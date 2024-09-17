'''
-*- coding: utf-8 -*-
Author: Parker Wise
Date: 2024-07-01
Description: Creates CSV of spectral cube properties from several files
Python Version: 3.11.9
'''
import astropy.io.fits as fits  # 6.1.0
from decimal import Decimal
import numpy as np  # 1.26.4

folderPath = "/home/pw/research/Cloud-C/co-data/"

fitsFilenames = ["A.Dust_Ridge_12C18O.cube.I.pbcor.fits",
                 "A.Dust_Ridge_13C16O.cube.I.pbcor.fits",
                 "A.Dust_Ridge_sci.spw29.cube.I.pbcor.fits",
                 "A.Dust_Ridge_sci.spw31.cube.I.pbcor.fits",
                 "B.Dust_Ridge_12C16O_1-0.cube.I.pbcor.fits",
                 "B.Dust_Ridge_12C17O.cube.I.pbcor.fits",
                 "B.Dust_Ridge_sci.spw29.cube.I.pbcor.fits",
                 "B.Dust_Ridge_sci.spw31.cube.I.pbcor.fits"]
# allows us to find header key given header value


def getHeaderList(header):
    try:
        headerKeys = np.array(list(header.keys()))
        headerValues = np.array(list(header.values()))
    except AttributeError:
        print("AttributeError: FITS file must have a header")
    return headerKeys, headerValues


def fetchFreqAxis(header):  # finds what axis freq is
    headerKeys, headerValues = getHeaderList(header)
    # it is assumed there is only one frequency axis
    freqIndex = np.where(  # or operators allows for all spectral types
            (headerValues == 'FREQ') |
            (headerValues == 'ENER') |
            (headerValues == 'WAVN') |
            (headerValues == 'VRAD') |
            (headerValues == 'WAVE') |
            (headerValues == 'VOPT') |
            (headerValues == 'ZOPT') |
            (headerValues == 'AWAV') |
            (headerValues == 'VELO') |
            (headerValues == 'BETA')
    )
    # returns false if there is no spectral axis
    freqKey = headerKeys[freqIndex]
    if len(freqKey) == 0:
        return False
    else:
        freqAxes = freqKey[0][-1]
    return freqAxes


def fetchSpatialAxes(header):
    headerKeys, headerValues = getHeaderList(header)
    spatialIndices = np.where(headerValues == 'deg')
    spatialKeys = headerKeys[spatialIndices]
    if len(spatialKeys) == 0:
        return False
    else:
        spatialAxes = [key[-1] for key in spatialKeys]
    return spatialAxes


def cube_properties(path):
    dataSet = path[0]
    path = folderPath+path
    image = fits.getdata(path)
    header = fits.getheader(path)
    spectralWindow = header['FILNAM04'][3:]
    xdim = header['NAXIS1']
    ydim = header['NAXIS2']
    spatialAxes = fetchSpatialAxes(header)
    if spatialAxes is False:
        pixelScale = "---"
    else:
        # pixel scale is multiplied from degrees to asec
        pixelScale = f"{Decimal(abs(header['CDELT1'])*3600):.2f}"
    bmaj = f"{Decimal(header['BMAJ']*3600):.2f}"
    bmin = f"{Decimal(header['BMIN']*3600):.2f}"
    bpa = f"{Decimal(header['BPA']):.2f}"
    freqAxis = fetchFreqAxis(header)
    if freqAxis is False:
        freqRange = "---"
        freqWidth = "---"
        channels = "---"
        channelWidth = "---"
    else:
        freqUnit = header['CUNIT'+freqAxis]
        # image is ordered in reverse: CTYPE4,CTYPE3...
        channels = np.shape(image)[-int(freqAxis)]
        # frequencies are converted to GHz instead of Hz
        startFreqVal = header['CRVAL'+freqAxis]*10**-9
        channelWidth = header['CDELT'+freqAxis]*10**-9
        endFreq = startFreqVal+channelWidth*channels
        freqWidth = (endFreq-startFreqVal)
        # Decimal formats values in scientific notation
        freqWidth = f"{Decimal(freqWidth):.2f}"
        startFreq = f"{Decimal(float(startFreqVal)):.2f}"
        endFreq = f"{Decimal(float(endFreq)):.2f}"
        freqRange = f"{startFreq} - {endFreq}"
        channelWidth = channelWidth * 10 ** 3
        channelWidth = f"{Decimal(channelWidth):.2f}"
    return startFreqVal, [f"{dataSet} & {spectralWindow} & {xdim} & {ydim}\\\\",
                          f"{dataSet} & {spectralWindow} & {pixelScale} & {bmaj} & {bmin} & {bpa} \\\\",
                          f"{dataSet} & {spectralWindow} & {channels} & {channelWidth} & {freqWidth} &  {freqRange}\\\\ \n"]

'''
csv should be sorted by ascending frequency
to avoid having to open all of our files twice, we save the freq and the
desired output to a dictionary that we then sort by frequency
'''
entryAndFreq = {}
for path in fitsFilenames:
    freq, entry = cube_properties(path)
    entryAndFreq[freq] = entry
sortedEntries = sorted(entryAndFreq.items())


def colString(numCols):
    colsString = "{|"
    for i in range(numCols):
        colsString += "c|"
    colsString += "}"
    return colsString


# Table 1
print(f"\\begin{{center}}\n\\begin{{tabular}}{colString(4)}")
print("Data Set& Spectral Window& X & Y \\\\ \n")
print("\\hline")
for entry in sortedEntries:
    print(entry[1][0])
print("\\end{tabular}\n\\end{center}")

# Table 2
print(f"\\begin{{center}}\n\\begin{{tabular}}{colString(6)}")
print("Data Set& Spectral Window &Pixel Scale(\"/px) & bmaj (\") & bmin (\") & bpa ($^{\\circ}$)\\\\ \n")
print("\\hline")
for entry in sortedEntries:
    print(entry[1][1])
print("\\end{tabular}\n\\end{center}")

# Table 3
print(f"\\begin{{center}}\n\\begin{{tabular}}{colString(6)}")
print("Data Set & spw & Channels & Channel Width (MHz) & Freq. Width (GHz) & Freq. Range (GHz) \\\\ \n")
print("\\hline")
for entry in sortedEntries:
    print(entry[1][2])
print("\\end{tabular}\n\\end{center}")
