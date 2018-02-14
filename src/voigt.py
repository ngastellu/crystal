#!/usr/bin/env python

from math import *
import numpy as np


def Voigt(xdata,xparams,gparam):
    """Takes in data and relevant parameters relevant to each point in the data
    and returns the Voigt profile evaluated at all data points."""

    F = np.zeros(len(xdata),dtype=float)
    
    for i in range(len(xdata)):
        amp = xparams[0][i]
        xc = xparams[1][i]
        w = xparams[2][i]
        
        gamma = (w/2.0)**2

        
        lorentzian = gamma/(np.power(xdata-xc,2)+gamma)
        gaussian = np.exp(-1.0*np.power(xdata-xc,2)/(2.0*(w/sqrt(2.0*log(2.0)))))

        F = np.add(amp*(gparam*lorentzian + (1-gparam)*gaussian),F,casting='unsafe') 

    return F
