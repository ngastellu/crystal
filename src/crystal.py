#!/usr/bin/env python


#idea: define lattice vectors and atomic positions within the unit cell in an input file that then gets fed to the calculation script
#IS IT NECESSARY? ==> maybe space group properties can be exploited to generate atomic positions within unit cell

#can get general reflections conditions from http://www.cryst.ehu.es/cryst/get_hkl.html for all 230 ptgroups

#so far this only takes into account general case i.e. special positions (reflection conditions associated to them) are ignored 

import numpy as np
from math import *
from atom import Atom

class Crystal:

    def __init__(self, buildFile):
        """Defines a crystal by its basis and its lattice vectors, which are themselves defined in
        the build file specified as an argument to this generator. This file also specifies whether or not
        asymmetries in electron density are to be taken into account or not."""

        self.Asymm = False #attributing determining whether or not to take asymmetries into account
        
        self.lattVs = np.zeros((3,3),dtype=float)
        basis = []

        with open(buildFile,'r') as fo:
            lines = fo.readlines()
        for line in lines:
            words = line.split()
            if words[0] == 'a':
                self.lattVs[0] = np.array(map(float,words[1:]))
            elif words[0] == 'b':
                self.lattVs[1] = np.array(map(float,words[1:]))
            elif words[0] == 'c':
                self.lattVs[2] = np.array(map(float,words[1:]))
            elif words[0] == 'atom':
                basis.append(Atom(words[1],words[2],map(int,words[3:])))
            elif words[0] == 'asymm':
                if words[1] == 'T':
                    self.Asymm = True
                elif words[1] == 'F':
                    self.Asymm = False
                else:
                    print "INPUT FILE ERROR: Only valid values for 'asymm' field are: 'T' ; 'F' .\nReturning 0.\n"
                    return 0
            else:
                print "INPUT FILE ERROR: Incorrect syntax\n" + line + "\nRefer to the README for correct input format.\nReturning 0\n"
                return 0

        #self.SpaceGroup = space_group #H-M symbol
        self.Basis = basis #array of Atom objects
    

    def getKlattice(self): #generates reciprocal lattice vecs using the real-spce latt vecs
        kVs = np.zeros((3,3),dtype=float)
        V = np.dot(self.lattVs[0],np.cross(self.lattVs[1],self.lattVs[2]))
        for i in range(2):
            kVs[i] = (2*pi/V)*np.cross(self.lattVs[i-2],self.lattVs[i-1])
        return kVs

    def getRConditions(self):
        """Returns an array strings of boolean expressions corresponding to the reflection conditions
        associated with the crystal's space group. 
        Assumes the script calling it labels the reciprocal lattice vector coefficients h,k,l."""
        
        space_group = self.SpaceGroup
        conditions = []

        if space_group == 'P21c': #one of many monoclinic space groups
            conditions.append('k==0 and l%2==0')
            conditions.append('h==0 and l==0 and k%2==0')
            conditions.append('h==0 and k==0 and l%2')

        elif space_group == 'I4/mcm': #for KCuF3
           conditions.append('(h+k+l)%2==0')
           conditions.append('l==0 and (h+k)%2==0')
           conditions.append('h==0 and k%2==0 and l%2==0')
           conditions.append('h==k and l%2==0')
           conditions.append('h==0 and k==0 and l%2==0')
           conditions.append('k==0 and l==0 and h%2==0')
        
        else:
           raise Exception("Space group %s not yet implemented" %(space_group))
    

    
