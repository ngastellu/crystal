#!/usr/bin/env python

import numpy as np
from math import *
from atom import Atom
from crystal import Crystal
from scipy.interpolate import interp1d, interp2d
from voigt import Voigt
import matplotlib.pyplot as plt

a = [8.270, 0.0, 0.0]
b = [0.0, 8.270, 0.0]
c = [0.0, 0.0, 7.841]

cop = Atom('Cu', '../data/Cu.xyz', [0,0,2,0,0])
pot = Atom('K','../data/K.xyz',[0,2])
fluo = Atom('F', '../data/F.xyz',[0,0,0,2])

system = Crystal([a,b,c],[pot,fluo,cop])
lattVs = system.lattVs #3x3 array containg the crystal's lattice vectors
kVs  = system.getKlattice() #3x3 array containing the reciprocal lattice vectors

for atm in system.Basis:
    atm.getAOparams() #loads fitting parameters necessary to calculation atomic form factors
    atm.getXYZ() #gets the coordinates of all atoms in unit cell
    atm.getSFa() #calculates atomic form factor for all atoms in basis

s_grid = np.array([0.001*k for k in range(1001)])
k_grid = np.array([4.0*pi*10.0*s for s in s_grid[:101]])
beta_grid = np.array([-pi + pi*k/50.0 for k in range(101)])
SKm = np.zeros((2,2197),dtype=complex) #matrix with each row containg h,k,l and |form factor|**2  
counter = 0
for h in [-i for i in range(-6,7)]:
    for k in [-i for i in range(-6,7)]:
        for l in [-i for i in range(-6,7)]:
            kVec = np.array(h*kVs[0]+k*kVs[1]+l*kVs[2])
            kNorm = np.linalg.norm(kVec)
            theta = np.arctan2(np.linalg.norm(np.cross(c,kVec)), np.dot(c,kVec))
            SFtmp = 0.0
            if kNorm <= k_grid[100] and kNorm!=0:
                print kNorm
                SKm[0][counter] = kNorm
                for atm in system.Basis:
                    ffs = np.zeros(len(atm.orbPop), dtype=float)
                    print np.shape(ffs)
                    for j in range(len(atm.XYZ)):
                        position = lattVs[0]*atm.XYZ[j][0] + lattVs[1]*atm.XYZ[j][1] + lattVs[2]*atm.XYZ[j][2]
                        if len(atm.orbPop) == 2: #dealing with only core and s AO
                            ffs = np.array([interp1d(k_grid,FF)(kNorm) for FF in atm.SFa])
                        elif len(atm.orbPop) == 4: #atom has core, s, and p(0-1) AOs.
                            ffs[:3] = np.array([interp1d(k_grid,FF)(kNorm) for FF in atm.SFa[:3]]) #need to do this case by case bc the nb of
                            ffs[3] = interp2d(k_grid, beta_grid, FF)(kNorm, theta) #symmetrical and asymmetrical AOs is different for each type of atom
                        else: #atom has core, s, and d(1-3) AOs
                            print np.shape(atm.SFa[2:])
                            ffs[:2] = np.array([interp1d(k_grid,FF)(kNorm) for FF in atm.SFa[:2]])
                            ffs[2:] = np.array([interp2d(k_grid,beta_grid,FF)(kNorm,theta) for FF in atm.SFa[2:]]).reshape(3)
                        ff = np.dot(atm.orbPop,ffs)
                        SFtmp += ff*np.exp(1j*np.dot(position,kVec))
            SKm[1][counter] = np.conj(SFtmp) * SFtmp
            print counter
            counter += 1

peak_widths = np.array([0.04 for j in range(2197)])
SKm = np.reshape(SKm,(2,2197))
print SKm
vparams = np.array([SKm[1],SKm[0],peak_widths])
print vparams
y = Voigt(4.0*pi*s_grid,vparams,0.5)

plt.plot(s_grid,y/np.amax(y), linewidth=1.0)
plt.xlim(0.1,0.83)
plt.ylim(0.0,1.1)
plt.show()

with open('../data/results_KCuF3rh_asymm.dat','w') as fo:
    for i in range(1001):
        fo.write('%s\t%s\n'%(str(s_grid[i]),str(y[i]/np.amax(y))))

