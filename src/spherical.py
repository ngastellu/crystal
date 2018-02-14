#!/usr/bin/env python

import numpy as np
from math import *
from atom import Atom
from crystal import Crystal

a = [0.0, 0.0, 5.720]
b = [4.53, 0.0, 0.0]
c = [0.0, 4.53, -2.86]

van = Atom('V', '../data/V.xyz', [0,0,0])
oxy1 = Atom('O','../data/O1.xyz',[0,0,0])
oxy2 = Atom('O', '../data/O2.xyz',[0,0,0])

system = Crystal([a,b,c],[van,oxy1,oxy2])
lattVs = system.LattVs #3x3 array containg the crystal's lattice vectors
kVs  = system.getKlattice() #3x3 array containing the reciprocal lattice vectors

for atm in system.Basis:
    atm.getAOparams() #loads fitting parameters necessary to calculation atomic form factors
    atm.getXYZ() #gets the coordinates of all atoms in unit cell
    atm.getSFs() #calculates atomic form factor for all atoms in basis

k_grid = 4.0*pi*van.SFs()[:,1]

for h in [-i for i in range(-6:7)]:
    for k in [-i for i in range(-6:7)]:
        for l in [-i for i in range(-6:7)]:
            Knorm = np.linalg.norm(h*kVs[0]+k*kVs[1]+l*kVs[2])








