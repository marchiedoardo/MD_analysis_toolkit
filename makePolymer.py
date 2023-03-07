#!/usr/bin/env python3

from math import sqrt,sin,cos
from scipy.stats import powerlaw,uniform,norm
import argparse
import numpy as np

def generate_ctcf_equid(N=1000, rho = 0.1):

   n_ctcf=0
   ctcf_list = []
   cnt = 1

   while n_ctcf<N*rho:

      ctcf_site = int(1/rho)*cnt 
      
      ctcf_type = 4
      ctcf_list.append([ctcf_site, ctcf_type])
      n_ctcf+=1
      cnt+=1

   return ctcf_list

def generate_ctcf_rand(N = 1000, rho = 0.1):
   
   n_ctcf = 0 
   ctcf_list = []
   
   while n_ctcf<N*rho:
      
      ctcf_site = int(uniform.rvs(1,N-1))
      if [ctcf_site,2] in ctcf_list or [ctcf_site,3] in ctcf_list or [ctcf_site,2] in ctcf_list:
         continue 
      ctcf_type = np.random.choice([2,3,4])
      ctcf_list.append([ctcf_site, ctcf_type])
      n_ctcf+=1
   return ctcf_list
      
def Parse():

   parser = argparse.ArgumentParser(description='Prepare simulations')
   parser.add_argument('-d', '--dataname', help='Lammps file data of the molecules') 
   parser.add_argument('-n', '--nmonomers', type=int, default=1000, help='Length of the polymer') 
   parser.add_argument('-b', '--blength', type=float, default=1, help='Monomer distance') 
   parser.add_argument('-m', '--mass', type=float, default=1000, help='Monomer mass') 
   parser.add_argument('-l', '--links', type=int, default=0, help='Number of random links') 
   parser.add_argument('-e', '--exponent', type=float, default=1, help='Exponent of the power-law for link connection') 
   parser.add_argument('-s', '--scale', type=float, default=-1, help='Scale of the power-law for link connection')
   parser.add_argument('-r', '--rho', type=float, default=0, help='Density of CTCF sites') 
   return parser.parse_args()

parms = Parse()
filename=parms.dataname
nMonomers=parms.nmonomers
bLength=parms.blength
mass=parms.mass
randomLinks=parms.links
exponent = parms.exponent
scale = parms.scale
rho = parms.rho

if rho != 0:
    ctcf_list = generate_ctcf_equid(nMonomers, rho)

k = bLength / sqrt( (cos(1)-1)**2 + sin(1)**2 + 1 )
side = k * nMonomers * 1.1

with open(filename, 'w') as f:

   print("LAMMPS Description", file=f)
   print("",file=f)
   print(nMonomers,"atoms",file=f)
   print(nMonomers-1+randomLinks,"bonds",file=f)
   print("",file=f)
   print("1 atom types",file=f)
   print("1 bond types",file=f)
   print("",file=f)
   print("%7.3f %7.3f xlo xhi" % (-side,side), file=f)
   print("%7.3f %7.3f ylo yhi" % (-side,side), file=f)
   print("%7.3f %7.3f zlo zhi" % (-side,side), file=f)
   print("",file=f)
   print("Masses",file=f)
   print("",file=f)
   print("1 ",mass,file=f)
   print("",file=f)
   
   print("Atoms",file=f)
   print("",file=f)
   for i in range(nMonomers):

      z = k * float(i) + k / 10
      x = k * ( cos(i) + 1.1 )
      y = k * ( sin(i) + 1.1 )
      
      if [i,2] in ctcf_list:
         print("%d %d %d\t %8.4f %8.4f %8.4f" % (i+1,1,2,x,y,z), file=f)
      elif [i,3] in ctcf_list:
         print("%d %d %d\t %8.4f %8.4f %8.4f" % (i+1,1,3,x,y,z), file=f)
      elif [i,4] in ctcf_list: 
         print("%d %d %d\t %8.4f %8.4f %8.4f" % (i+1,1,4,x,y,z), file=f)
      else: 
         print("%d %d %d\t %8.4f %8.4f %8.4f" % (i+1,1,1,x,y,z), file=f)

   print("",file=f)
   print("Bonds",file=f)
   print("",file=f)
   for i in range(nMonomers-1):

      print(i+1,1,i+1,i+2,file=f)

   if randomLinks>0:

     if scale<0:
       scale = nMonomers

     i = 0
     with open('tmplist.dat', 'w') as g:

       while i<randomLinks:
          
         #dist = int(powerlaw.rvs(exponent, scale=scale))+1
         dist = int(norm.rvs(loc = 200, scale = 10))
         start = int(uniform.rvs(1,nMonomers-1))
         #end = int(uniform.rvs(1,nMonomers-1))
         if (start+dist<nMonomers and start!=(start+dist)):
            print(i+nMonomers,1,start,start+dist,file=f)
            #print(start,start+dist,file=g)
            i=i+1
         #print(i+nMonomers, 1, start, end , file=f)
         #i=i+1

#with open('ctcf_sites_rnd.data', 'w') as h:
 #  for i in range(len(ctcf_list)):
  #     print(ctcf_list[i][0], ctcf_list[i][1], file=h)

