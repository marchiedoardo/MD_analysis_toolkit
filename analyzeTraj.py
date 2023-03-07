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

def Distance2(a, b):
   return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2
 
def Plot(timeArray, dMsd, figName, fit, dMsdCom=np.array([]), fitCom=np.array([])):
   plt.loglog(timeArray, dMsd,"b", label="bead")
   plt.loglog(timeArray, fit, "b--")
   if dMsdCom.size != 0:
     plt.loglog(timeArray, dMsdCom, "r", label="com")
   if fitCom.size != 0:
     plt.loglog(timeArray, fitCom, "r--")
   plt.legend(loc="upper left")
   plt.xlabel("t")
   plt.ylabel("MSD")
   if figName=='screen':
      plt.show()
   else:
      plt.savefig(figName)

#computing mean square displ.
def CalculateDmsd(pos):

   nFrames = len(pos)
   dSum = np.zeros(nFrames)
   count = np.zeros(nFrames)
   dMsd = np.zeros(nFrames)

   for i in tqdm(range(nFrames)):
      for j in range(i,nFrames):
        d = Distance2(pos[i], pos[j])
        di = j - i
        dSum[di] = dSum[di] + d
        count[di] = count[di] + 1

   for di in range(nFrames):
      dMsd[di] = dSum[di] / count[di]

   return dMsd

def Fit(timeArray, data, fitLimit=-1):

   x = np.log(timeArray[1:])
   y = np.log(data[1:])
           
   if fitLimit>0:
      x = x[:fitLimit]
      y = y[:fitLimit]

   lr = linregress(x, y)
   #print(lr)
   alpha = lr.slope
   diffCoeff = np.exp(lr.intercept)
   pValue = lr.pvalue
   fit = diffCoeff * np.power( timeArray, alpha )

   print("alpha=",alpha)
   print("D=",diffCoeff)
   print("p-value=",pValue)

   return fit

def Parse():

   parser = argparse.ArgumentParser(description='Analyze MDE data')
   parser.add_argument('-f', '--trajname', help='File dcd of the trajectory') 
   parser.add_argument('-d', '--dataname', help='Lammps file data of the molecules') 
   parser.add_argument('-n', '--natom', type=int, help='What atom to follow') 
   parser.add_argument('-s', '--skip', type=int, default = 0, help='Skip first frames')
   parser.add_argument('-e', '--deltat', type=float, default = 1, help='Number of timesteps between frames') 
   parser.add_argument('-l', '--fitlimit', type=int, default=-1, help='Fit only up to a limit') 
   parser.add_argument('-t', '--truncate', type=int, default=-1, help='truncate the traj at given element') 
   parser.add_argument('-b', '--beadname', default='none', help='Filename to print bead data') 
   return parser.parse_args()
   
##############################################################

parms = Parse()

nAtom = parms.natom
fileName = parms.trajname
dataName = parms.dataname
fitLimit = parms.fitlimit
deltat = int(parms.deltat)
skip = int(parms.skip)
truncate = parms.truncate

beadName = parms.beadname

# Load configuration and trajectory with MDanalysis 
u = md.Universe(dataName, fileName, format="LAMMPS")
print("File=",fileName)
dt = u.trajectory.dt * deltat 
tmax = u.trajectory.totaltime * deltat 

print("Tmax=",tmax," (npoints=",tmax/dt,")")
print("dt=",dt)

#Choosing atom to follow 
selAtom = u.select_atoms("index "+str(nAtom))
#all the atoms?
selAll = u.residues.atoms

#get positions of the selected atom and of the center of mass
previousCoord = selAtom.positions
previousCom = selAll.center_of_mass()

timeArray = np.arange( 0, tmax + dt, dt )

# Read from trajectory

posTraj = []
posCom = []

for ts in u.trajectory:

   pos = selAtom.positions[0]
   posTraj.append(pos)
   posCom.append(selAll.center_of_mass())

if truncate>0:
   posTraj = posTraj[:truncate]
   posCom = posCom[:truncate]
   timeArray = timeArray[:truncate]

if stride>1:
   posTraj = posTraj[::stride]
   posCom = posCom[::stride]
   timeArray = timeArray[::stride]

# Calculate stuff

print("Calculate MSD of a bead")

posTraj = posTraj[skip:]
PosCom = posCom[skip:]

dMsd = CalculateDmsd(posTraj)

print()
print("Bead #",nAtom)

if beadName != 'none':
    with open(beadName, 'w') as f:
       for i, d in enumerate(dMsd):
          print( timeArray[i], d, file=f )

