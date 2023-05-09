#!/usr/bin/python

#/usr/bin/env python3.9

import MDAnalysis as md
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib import cm
from tqdm import tqdm
import argparse

def radgyr(atoms, masses, total_mass=None):
    coordinates = atoms.positions
    center_of_mass = atoms.center_of_mass()
    
    #get squared distances
    ri_sq = (coordinates-center_of_mass)**2
    #sum positions
    sq = np.sum(ri_sq, axis=1)
    #weight positions
    rog_sq = np.sum(masses*sq, axis=1)/total_mass

    return np.sqrt(rog_sq)

def Parse():

   parser = argparse.ArgumentParser(description='Analyze MDE data')
   parser.add_argument('-f', '--trajname', help='File dcd of the trajectory') 
   parser.add_argument('-d', '--dataname', help='Lammps file data of the molecules') 
   parser.add_argument('-b', '--rogname', default='none', help='Filename to print rog data')
   parser.add_argument('-e', '--deltat', type=float, help='Time difference between conf in the traj file')
   return parser.parse_args()
   
  
##############################################################

#Inputs
parms = Parse()

deltat = parms.deltat
rogName = parms.rogname
fileName = parms.trajname
dataName = parms.dataname

# Load configuration and trajectory with MDanalysis 
u = md.Universe(dataName, fileName, format="LAMMPS")
print("File=",fileName)

dt = u.trajectory.dt * deltat
tmax = u.trajectory.totaltime*deltat

timeArray = np.arange( 0, tmax + dt, dt )

#Select all atoms
selAll = u.residues.atoms

rogTraj = []

#Compute rog
for ts in u.trajectory:
   rog = selAll.radius_of_gyration()
   rogTraj.append(rog)

#Write rog data file
if rogName != 'none':
    with open(rogName, 'w') as f:
       for i, d in enumerate(rogTraj):
          print( timeArray[i], d, file=f )

