#!/usr/bin/env python

#*****not yet implemented for H and He***** 
import numpy as np
from math import *

class Atom:

    def __init__(self,name,orb_pop): #define atom based on its name and a list of its vance AOs populations
                                   #which is loaded as a list if ints, the ordering of valence oribitals 
                                   #in this list is the same as the one listed in 'ffs_sorted.dats'
        self.name = name
        self.orbPop = orb_pop
    
    def getAOparams(self): #obtain Gaussian fit parameters for aspherical electron scattering of the atom's different
                           #valence orbitals
        with open('../data/ffs_sorted.dat','r') as fo:
            self.AOparams = []
            flines = fo.readlines()
            for line in flines:
                if line.split()[0]==self.name:
                    self.AOparams.append(map(float,line.split()[3:])) #creates an array of the form[[a params for spherical FF],[b params for spherical FF],
                                                                    #[a's for core],[b's for core],...[a's for last AO],[b's for the last AO]]
            self.AOparams = np.array(self.AOparams)    
    
    def getSFs(self): #calculates the atomic spherical structure
        s_grid = np.array([0.01*k for k in range(101)])
        self.SFs = np.zeros(101, dtype=float)
       
        for k in range(8):
            self.SFs = np.add(self.AOparams[0][k]*np.exp(-1.0*self.AOparams[1][k]*np.square(s_grid)), self.SFs)


   
    def getSFa(self): #calculates the aspherical structure factor for atom in question
         s_grid = np.array([0.01*k for k in range(101)])
        
         dim = len(self.orbPop)
        
        
         if dim == 1: #dealing with only s AO
             self.SFa = np.zeros(101,dtype=float)
             for k in range(8):
                 self.SFa = np.add(self.AOparams[-2][k]*np.exp(-self.AOparams[-1][k]*np.square(s_grid)), self.SFa) #
        
        
         elif dim == 3: #dealing with s, p0, and p1 AOs
             FF_s = np.zeros(101,dtype=float)
             FF_p0 = np.zeros(101,dtype=float)
             ff_p1 = np.zeros(101,dtype=float)
             FF_p1 = np.zeros((101,101),dtype=float)
            
             for k in range(8):
                 FF_s = np.add(self.AOparams[4][k]*np.exp(-1.0*self.AOparams[5][k]*np.square(s_grid)), FF_s)
                 FF_p0 = np.add(self.AOparams[6][k]*np.exp(-1.0*self.AOparams[7][k]*np.square(s_grid)), FF_p0) 
                 ff_p1 = np.add(self.AOparams[8][k]*np.exp(-1.0*self.AOparams[9][k]*np.square(s_grid)), ff_p1) 
             beta_grid = np.array([-pi+k*pi/50.0 for k in range(101)])
             a_p = (3.0/2.0)*np.square(np.sin(beta_grid))
             b_p = np.square(np.cos(beta_grid)) - (0.5)*np.square(np.sin(beta_grid))
             for k in range(101):
                 FF_p1[k] = np.add(a_p*FF_p0[k],b_p*ff_p1[k])
             
             self.SFa = [FF_s,FF_p0,FF_p1]
            


         else: #dealing with s, d1, d2, d3 AOs
             FF_s = np.zeros(101,dtype=float)
             ff_d1 = np.zeros(101,dtype=float)
             ff_d2 = np.zeros(101,dtype=float)
             ff_d2 = np.zeros(101,dtype=float)
             ff_d3 = np.zeros(101,dtype=float)
             FF_d1 = np.zeros((101,101),dtype=float)
             FF_d2 = np.zeros((101,101),dtype=float)
             FF_d3 = np.zeros((101,101),dtype=float)
            
             for k in range(8):
                 FF_s = np.add(self.AOparams[4][k]*np.exp(-self.AOparams[5][k]*np.square(s_grid)), FF_s)
                 ff_d1 = np.add(self.AOparams[6][k]*np.exp(-self.AOparams[7][k]*np.square(s_grid)), ff_d1) 
                 ff_d2 = np.add(self.AOparams[8][k]*np.exp(-self.AOparams[9][k]*np.square(s_grid)), ff_d2)
                 ff_d3 = np.add(self.AOparams[10][k]*np.exp(-self.AOparams[11][k]*np.square(s_grid)), ff_d3)

             beta_grid = np.array([-pi+k*(pi/50.0) for k in range(101)])
            
             a_d1 = (1.0/4.0) + 1.50*np.square(np.cos(beta_grid)) + (9.0/4.0)*np.power(np.cos(beta_grid),4)
             b_d1 = 3.0*(np.square(np.cos(beta_grid)) - np.power(np.cos(beta_grid),4))
             c_d1 = (3.0/4.0) - 1.50*np.square(np.cos(beta_grid)) + 0.75*np.power(np.cos(beta_grid),4)
            
             a_d2 = 1.50*(np.square(np.cos(beta_grid)) - np.power(np.cos(beta_grid),4))
             b_d2 = 0.50 - 1.50*np.square(np.cos(beta_grid)) + 2.0*np.power(np.cos(beta_grid),4)
             c_d2 = 0.50*(1 - np.power(np.cos(beta_grid),4))
            
             a_d3 = (3.0/8.0) - 0.75*np.square(np.cos(beta_grid)) + (3.0/8.0)*np.power(np.cos(beta_grid),4)
             b_d3 = 0.50*(1.0 - np.power(np.cos(beta_grid),4))
             c_d3 = (1.0/8.0) + 0.75*np.square(np.cos(beta_grid)) + (1.0/8.0)*np.power(np.cos(beta_grid),4)

             for k in range(101):
                 FF_d1[k] = a_d1*ff_d1[k] + b_d1*ff_d2[k] + c_d1*ff_d3[k]
                 FF_d2[k] = a_d2*ff_d1[k] + b_d2*ff_d2[k] + c_d2*ff_d3[k]
                 FF_d3[k] = a_d3*ff_d1[k] + b_d3*ff_d2[k] + c_d3*ff_d3[k]
             
             self.SFa = [FF_s,FF_d1,FF_d2]
        
                
                                   

                
