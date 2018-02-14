#!/usr/bin/env python

import numpy as np
from math import *
from atom import Atom
from crystal import Crystal
from scipy.interpolate import interp1d
from voigt import Voigt
import matplotlib.pyplot as plt

a = [0.0, 0.0, 5.720]
b = [4.53, 0.0, 0.0]
c = [0.0, 4.53, -2.86]

van = Atom('V', '../data/V.xyz', [0,0,0])
oxy1 = Atom('O','../data/O1.xyz',[0,0,0])
oxy2 = Atom('O', '../data/O2.xyz',[0,0,0])

system = Crystal([a,b,c],[van,oxy1,oxy2])
lattVs = system.lattVs #3x3 array containg the crystal's lattice vectors
kVs  = system.getKlattice() #3x3 array containing the reciprocal lattice vectors

for atm in system.Basis:
    atm.getAOparams() #loads fitting parameters necessary to calculation atomic form factors
    atm.getXYZ() #gets the coordinates of all atoms in unit cell
    atm.getSFs() #calculates atomic form factor for all atoms in basis

s_grid = np.array([0.001*k for k in range(1001)])
k_grid = np.array([4.0*pi*10.0*s for s in s_grid[:101]])
SKm = np.zeros((216,2),dtype=float) #matrix with each row containg h,k,l and |form factor|**2  
counter = 0
for h in [-i for i in range(-6,7)]:
    for k in [-i for i in range(-6,7)]:
        for l in [-i for i in range(-6,7)]:
            kVec = np.array(h*kVs[0]+k*kVs[1]+l*kVs[2])
            kNorm = np.linalg.norm(kVec)
            SFtmp = 0.0
            if kNorm <= k_grid[100] and kNorm!=0:
                SKm[counter][0] = kNorm
                for atm in system.Basis:
                    for j in range(len(atm.XYZ)):
                        position = lattVs[0]*atm.XYZ[j][0] + lattVs[1]*atm.XYZ[j][1] + lattVs[2]*atm.XYZ[j][2]
                        ff = interp1d(k_grid,atm.SFs)(kNorm)
                        SFtmp += ff*np.exp(1j*np.dot(position,kVec))
            
            SKm[counter][1] = np.conj(SFtmp) * SFtmp
            counter += 1

peak_widths = np.array([0.12 for j in range(216)])

vparams = np.array([SKm[:][1],SKm[:][0],peak_widths])
y = Voigt(4.0*pi*s_grid,vparams,0.5)

plt.plot(s_grid,y/np.amax(y))
plt.show()

                

                            
