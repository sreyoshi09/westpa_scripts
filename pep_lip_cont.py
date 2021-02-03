#!/usr/bin/env python

import loos
import loos.pyloos
import sys
import numpy
import copy
import os
import h5py
import argparse

parser = argparse.ArgumentParser()
#### usually used when cholesterol is not present in the dataset
### Calculates the fengycin-lipid contacts at each frame for each trajectory in the WE simulations. Outputs all this as a hdf file

## System details from command line
parser.add_argument('dir',help='Path of the system', type=str)
parser.add_argument('model', help='PSF of the system',type=str)
parser.add_argument('feng', help ='select the atoms for fengycin',type=str)
parser.add_argument('lipid', help ='select the lipid atoms',type=str)
parser.add_argument('first',help = 'start analysis from this iteration', type=int)
parser.add_argument('last',help = 'end analysis at this iteration', type=int)

args = parser.parse_args()


##Select only the fengycin tail and lipid tails

system = args.dir
source = loos.createSystem(args.model)
feng = loos.selectAtoms(source,args.feng)

lipid = loos.selectAtoms(source, args.lipid)
individual_lipid = lipid.splitByMolecule()
## Where the trajectory files are stored and westpa config file is stored
root_dir=system + "/"+"traj_segs"
westh5 = system + "/" + "west.h5"
iter = []

### Name of the output file
outputh5file = system + "/"+ "feng_lipid.h5"
iter_sorted = []
for x in os.listdir(root_dir):    
    if (int(x) >= args.first and int(x)<=args.last):
        iter.append(x)
iter_sorted = sorted(iter,key=int)
iterations = numpy.array([])
### Create an object of the hdf file
f=h5py.File(outputh5file,"w")
g= f.create_group("iterations")
west =h5py.File(westh5,"r")
### This function calculates the fengycin-lipid contacts as a function of distance between fengycin and each lipid
def lipid_feng(source,feng,individual_lipid,box):
    upper_leaflet = loos.AtomicGroup()
#    feng.reimage()
    individual_feng = feng.splitByMolecule()

    for mol in individual_lipid:
        phosphates = loos.selectAtoms(mol, 'name == "PO4"')
        if (phosphates[0].coords()).z() > 0.0:
            upper_leaflet.append(mol)
    CAB_fl=0
    upper_leaflet_C = loos.selectAtoms(upper_leaflet, "name =~ '^[CD][0-9][AB]?'")
### Only the upper leaflet lipids are considered because fengycins are also in the upper leaflet.
    indiv_u_l = upper_leaflet.splitByMolecule()
    for i in indiv_u_l:
        for j in individual_feng:
            cent_1=i.centroid()
            cent_2=j.centroid()
            distance_between1 = cent_1.distance(cent_2,box)
            Sij_fl = 1/(1+(distance_between1/10.0)**6)
            CAB_fl += Sij_fl
                    
    return CAB_fl



#### h5py architecture specific code. Going through each WESTPA iterations, then each bin of the progress coordinate and finally applying the feng-lipid function on each frame of the trajectory present in that bin.
for i in iter_sorted:
    t = int(i)
    iter_dir = root_dir + '/' + i
    seg=[]
    data = numpy.empty([len(os.listdir(iter_dir)),6])  
    for j in os.listdir(iter_dir):
        seg.append(j)
    seg_dir=sorted(seg,key=int)
    for c in seg_dir:
        l=int(c)
        seg = root_dir + '/' + i + '/' + c+ '/'+'seg.dcd'
        traj = loos.pyloos.Trajectory(seg,source)
        k=0
        for frame in traj:
           box = source.periodicBox()
           data[l][k]=lipid_feng(source,feng,individual_lipid,box)
           print i,'\t',c,'\t',data[l][k]
           k = k+1
##  Specific to the nomencalture in west.cfg file          
    name_ds='iter'+'_'+'00'+i
    name_subgroup = 'iterations'+'/'+name_ds+'/'+'seg_index' 
    pcoord_name = 'iterations'+'/'+name_ds+'/'+'pcoord'
    pcoords_data = west[pcoord_name]
    pcoord = numpy.array(pcoords_data)
    final = numpy.empty([pcoord.shape[0],pcoord.shape[1],2])

    for u in range(0, pcoord.shape[0]):
        for v in range(0,pcoord.shape[1]):
            final[u][v][0]=pcoord[u][v][0]
            final[u][v][1]=data[u][v]                           
## Specified variable name feng-lipid is called in the add_on.py file.                
    s=g.create_group(name_ds)
    s.attrs['n_iter']= t
    s.create_dataset("feng-lipid",shape=(len(os.listdir(iter_dir)),6,2),data=final)
west.close()
f.close()
