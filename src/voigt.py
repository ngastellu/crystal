#!/usr/bin/env python

from math import *
import numpy as np


def Voigt(xdata,xparams,gparam):
    """Takes in data and relevant parameters relevant to each point in the data
    and returns the Voigt profile evaluated at all data points."""

    F = np.zeros(len(xdata),dtype=float)
    
    for i in range(len(xdata)):
        amp = xparams[i][0]
        x0 = xparams[i][1]
        w = xparams[i][2]

        F += gparam*amp*(((w/2.0)**2)*(np.power((xdata-xc),-2)+(w/2.0)**2)) + (1-gparam)*amp*np.exp(np.power((-1.0*(xdata-xc),2))/(2.0*(w/(2.0*sqrt(2.0*log(2.0))))**2))
    
    return F
