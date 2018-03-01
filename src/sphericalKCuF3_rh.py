#!/usr/bin/env python

import numpy as np
from math import *
from atom import Atom
from crystal import Crystal
from scipy.interpolate import interp1d
from voigt import Voigt
import matplotlib.pyplot as plt

a = [8.270, 0.0, 0.0]
b = [0.0, 8.270, 0.0]
c = [0.0, 0.0, 7.841]

pot = Atom('K', '../data/K.xyz',[0,0,0])
cop = Atom('Cu', '../data/Cu.xyz', [0,0,0])
flu = Atom('F','../data/F.xyz',[0,0,0])

system = Crystal([a,b,c],[pot,cop,flu])
lattVs = system.lattVs #3x3 array containg the crystal's lattice vectors
kVs  = system.getKlattice() #3x3 array containing the reciprocal lattice vectors

for atm in system.Basis:
    atm.getAOparams() #loads fitting parameters necessary to calculation atomic form factors
    atm.getXYZ() #gets the coordinates of all atoms in unit cell
    atm.getSFs() #calculates atomic form factor for all atoms in basis

s_grid = np.array([0.001*k for k in range(1001)])
k_grid = np.array([4.0*pi*10.0*s for s in s_grid[:101]])
SKm = np.zeros((2,2197),dtype=complex) #matrix with each row containg h,k,l and |form factor|**2  
counter = 0
for h in [-i for i in range(-6,7)]:
    for k in [-i for i in range(-6,7)]:
        for l in [-i for i in range(-6,7)]:
            kVec = np.array(h*kVs[0]+k*kVs[1]+l*kVs[2])
            kNorm = np.linalg.norm(kVec)
            SFtmp = 0.0
            if kNorm <= k_grid[100] and kNorm!=0:
                print kNorm
                SKm[0][counter] = kNorm
                #print 'ye'
                for atm in system.Basis:
                    for j in range(len(atm.XYZ)):
                        position = lattVs[0]*atm.XYZ[j][0] + lattVs[1]*atm.XYZ[j][1] + lattVs[2]*atm.XYZ[j][2]
                        ff = interp1d(k_grid,atm.SFs)(kNorm)
                        SFtmp += ff*np.exp(1j*np.dot(position,kVec))
            
            SKm[1][counter] = np.conj(SFtmp) * SFtmp
            print counter
            counter += 1

peak_widths = np.array([0.12 for j in range(2197)])
SKm = np.reshape(SKm,(2,2197))
print SKm
vparams = np.array([SKm[1],SKm[0],peak_widths])
print vparams
y = Voigt(4.0*pi*s_grid,vparams,0.5)

plt.plot(s_grid,y/np.amax(y))
plt.show()

with open('../data/results_VO2m.dat','w') as fo:
    for i in range(1001):
        fo.write('%s\t%s\n'%(str(s_grid[i]),str(y[i]/np.amax(y))))

