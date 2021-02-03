#!/usr/bin/env python

import loos
import loos.pyloos
import sys
import numpy
import copy
import os
import h5py
import argparse


#### This is a script that calculates number of aggregtes of fengycin bound to membrane and stores it in an hdf file
### This will be used by w_pdist in WESTPA to calculate the probability distribution
parser = argparse.ArgumentParser()

parser.add_argument('dir',help='Path of the system', type=str)
parser.add_argument('model', help='PSF of the system', type=str)
parser.add_argument('feng', help='select atoms for fengycin')
parser.add_argument('first',help = 'start analysis from this iteration', type=int)
parser.add_argument('last',help = 'end analysis at this iteration', type=int)
args = parser.parse_args()

### function calculates the no. of aggregates
def no_of_aggs(molecules,box): 
###Makes two lists molecules copy with each fengycin a member of this list and clusters where we keep adding the fengycins which are part of an aggregate
    molecules_copy = list(molecules)
    clusters=[molecules_copy.pop(0)]
#   print "function starts",type(clusters)
    ###p is the counter
    p = 0
    #### This list goes on till all molecules are assigned to a particular aggregate in the list clusters
    while(len(molecules_copy)>0):
        non_contacts = []
        flags=0
        for mol in molecules_copy:
            if(mol.contactWith(6.5,clusters[p],box,20)==True):
#               print "Inside loop",type(clusters)
#               print "element of cluster",type(clusters[p])
                ####Part of one cluster so appended to that
                clusters[p].append(mol)
                flag=1
            else:
                #### Not a part of clusters are assigned to single lists non_contacts
                non_contacts.append(mol)
        ###AFter going through all the molecules which are not part of any cluster the non_contacts one become new molecules_copy
        molecules_copy=non_contacts
        if(flags!=1):
            #### If a molecule is found part of the aggregate the first fengycin molecule is moved to clusters from molecules and the cycle continues
            if(len(molecules_copy)>0):
                p=p+1
                clusters.append(molecules_copy.pop(0))
    print  len(clusters)
    return len(clusters)


### System specific requirements
system = args.dir
psffile =args.model
source = loos.createSystem(psffile)
fengs = loos.selectAtoms(source,args.feng)
molecules= fengs.splitByMolecule()

### Go through the file system again and gain access to west.h5
root_dir=system + "/"+"traj_segs"
westh5 = system + "/" + "west.h5"
iter = []
outputh5file = system + "/"+ "agg.h5"
iter_sorted = []
for x in os.listdir(root_dir):
    if (int(x) >= args.first and int(x)<=args.last):
        iter.append(x)
 
iter_sorted = sorted(iter,key=int)
iterations = numpy.array([])
f=h5py.File(outputh5file,"w")
g= f.create_group("iterations")
westh5 = system + "/" + "west.h5"
west = h5py.File(westh5,"r")
#calculat the number of aggregates and put them in a numpy array

for i in iter_sorted:
    
    iter_dir = root_dir + '/' + i
    t = int(i)
    seg=[]
    
    data = numpy.zeros([len(os.listdir(iter_dir)),6])  
    for j in os.listdir(iter_dir):
        seg.append(j)
    seg_dir=sorted(seg,key=int)
    for c in seg_dir:
         l=int(c)
         seg = root_dir + '/' + i + '/' + c+ '/'+'seg.dcd'
         traj = loos.pyloos.Trajectory(seg,source)
         k=0
         for frame in traj:
            molecules= fengs.splitByMolecule()  
            box = source.periodicBox()
            print i,"\t",c
            ### Call the no. of aggs function
            data[l][k] = no_of_aggs(molecules,box) 
            print  data[l][k]
            k=k+1
#retrieve the pcoord values
#### Store both the progress coordinate from west.h5 and the no. of aggs in a single h5 file
    name_ds='iter'+'_'+'00'+i
    pcoord_name = 'iterations'+'/'+name_ds+'/'+'pcoord'
    pcoords_data = west[pcoord_name]
    pcoord = numpy.array(pcoords_data)
    final = numpy.empty([pcoord.shape[0],pcoord.shape[1],2])

    for u in range(0, pcoord.shape[0]):
        for v in range(0,pcoord.shape[1]):
            final[u][v][0]=pcoord[u][v][0]
            final[u][v][1]=data[u][v]
            print final[u][v][1]

    name_subgroup = 'iterations'+'/'+name_ds+'/'+'seg_index'     
    s=g.create_group(name_ds)
    s.attrs['n_iter']= t
    s.create_dataset("aggs",shape=(len(os.listdir(iter_dir)),6,2),data=final)
f.close()
