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
### also generates weighted average distribution of mers


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
    ###p is the counter
    p = 0
    #### This list goes on till all molecules are assigned to a particular aggregate in the list clusters
    while(len(molecules_copy)>0):
        non_contacts = []
        flags=0
        for mol in molecules_copy:
            if(mol.contactWith(7,clusters[p],box,1200)==True):
                ####Part of one cluster so appended to that 
                clusters[p].append(mol)
                flag=1
            else:
                #### Not a part of clusters are assigned to single lists non_contacts
                non_contacts.append((molecules_copy.pop(0)))
        ###AFter going through all the molecules which are not part of any cluster the non_contacts one become new molecules_copy
        molecules_copy=non_contacts
        if(flags!=1):
            #### If a molecule is found part of the aggregate the first fengycin molecule is moved to clusters from molecules and the cycle continues
            if(len(molecules_copy)>0):
                p=p+1
                clusters.append(molecules_copy.pop(0))
    print  len(clusters)
    return clusters

def no_of_mers(clusters):
   no_of_mers = numpy.zeros((15,1))
   for i in range(len(clusters)):
       mers=clusters[i].splitByMolecule()
        print mers
       for j in range(0,15):
           if len(mers)==j:
               no_of_mers[j]+=1

   print "no. of mers","\t",no_of_mers
   print  len(clusters)
   return no_of_mers


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
outputh5file = system + "/"+ "agg_1d.h5"


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
#mers_avg =  numpy.zeros((15))
#mers_frame = []
#traj_wt = []
for i in iter_sorted:
    
    iter_dir = root_dir + '/' + i
    t = int(i)
    seg=[]
    name_ds='iter'+'_'+'00'+i    
    sg_name = "iterations/iter_00" + i + "/seg_index"
    sg_data = west[sg_name]
    seg_index = numpy.array(sg_data)
    data = numpy.zeros([len(os.listdir(iter_dir)),6])  
    for j in os.listdir(iter_dir):
        seg.append(j)
    seg_dir=sorted(seg,key=int)
    for c in seg_dir:     
         l=int(c)
         wt = seg_index[l][0] 
         seg = root_dir + '/' + i + '/' + c+ '/'+'seg.dcd' 
         traj = loos.pyloos.Trajectory(seg,source) 
         k=0 
         for frame in traj: 
            molecules = fengs.splitByMolecule() 
#            print molecules
            box = source.periodicBox() 
            print i,"\t",c ### Call the no. of aggs function
            
            clusters = no_of_aggs(molecules,box) 
            data[l][k] = len(clusters)
            # print data[l][k]
            
#           # if len(clusters) == 3:
#           traj_wt.append(wt) 
#           mers_frame.append(no_of_mers(clusters))
#           # print  data[l][k]
            k=k+1
#
    name_subgroup = 'iterations'+'/'+name_ds+'/'+'seg_index'
    s=g.create_group(name_ds) 
    s.attrs['n_iter']= t
    s.create_dataset("aggs_1d",shape=(len(os.listdir(iter_dir)),6,1),data=data)

f.close()

#um_traj_wt = sum(traj_wt)
#for stuff in range(len(mers_frame)):
#    for value in range(len(mers_frame[stuff])):
#        mers_avg[value] += mers_frame[stuff][value]*traj_wt[stuff]*(value+1)/(sum_traj_wt*15)
#
#for mer in range(len(mers_avg)):
#    print mers_avg[mer]
