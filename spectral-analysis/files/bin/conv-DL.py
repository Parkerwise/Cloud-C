#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Compute conventional Drude-Lorentz model
#
#
#  Who                  When            What
#
#  T. Moeller           2024-02-11      initial version
#
#******************************************************************************


#------------------------------------------------------------------------------
# import python packages
import os
import numpy as np
#------------------------------------------------------------------------------


##-----------------------------------------------------------------------------
##
## start main program
##
if __name__ == '__main__':


    ## import input file "in.txt"
    InFileName = "in.txt"
    InFile = open(InFileName, 'r')
    InFileContent = InFile.readlines()
    InFile.close()
    w0 = []
    wp = []
    G = []
    LineCounter = (-1)
    for LineID, line in enumerate(InFileContent):
        StrippedLine = line.strip()
        if (StrippedLine != ""):
            LineCounter += 1


            ## read number of frequency points
            if (LineCounter == 0):
                NumberXValues = int(StrippedLine)


            ## value of epsilon infinity
            elif (LineCounter == 1):
                Einf = float(StrippedLine)


            ## read number of oscillators
            elif (LineCounter == 2):
                number_osc = int(StrippedLine)

            else:
                SplittedLine = StrippedLine.split()
                w0.append(float(SplittedLine[0]))
                wp.append(float(SplittedLine[1]))
                G.append(float(SplittedLine[2]))


    ## convert to numpy array
    w0 = np.asarray(w0)
    wp = np.asarray(wp)
    G = np.asarray(G)


    ## get input file "DataIn.dat"
    DataInFileName = "DataIn.dat"
    DataInFileContent = np.loadtxt(DataInFileName)


    ## compute reflectance
    Reflectance = []
    for i in range(NumberXValues):
        w = DataInFileContent[i]
        eps = complex(0.0, 0.0)
        for j in range(1, number_osc + 1):
            if w == 0.0 and w0[j - 1] == 0.0:
                eps += (wp[j - 1]**2) / (w0[j - 1]**2 - (w + 0.00001)**2 - complex(0.0, 1.0) * G[j - 1] * (w + 0.00001))
            else:
                eps += (wp[j - 1]**2) / (w0[j - 1]**2 - w**2 - complex(0.0, 1.0) * G[j - 1] * w)
        eps = complex(Einf, 0.0) + eps
        refl = abs((1.0 - np.sqrt(eps)) / (1.0 + np.sqrt(eps)))**2
        Reflectance.append(refl)
    Reflectance = np.asarray(Reflectance)


    ## export reflectance
    OutputFileName = "FitFunctionValues.dat"
    Reflectance = np.savetxt(OutputFileName, Reflectance)


## finished!
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
