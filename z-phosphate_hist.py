#!/usr/bin/env python

import loos
import loos.pyloos
import sys
import numpy
import copy
import math
import os
import h5py
import argparse

####This script computes the difference in z-height of phosphates in the upper and lower leaflet for membranes with and without fengycins. Then it histograms the height of the phosphates w.r.t the distance between the phosphates and the fengycins.
###This is an all to all mapping. In this case we have incorporated the weights of each trajectory instead of using add_on.py.  

parser = argparse.ArgumentParser()

## All the necessary inputs for the script.
parser.add_argument('dir',help='Path of the system', type=str)
parser.add_argument('model', help='PSF of the system',type=str)
parser.add_argument('feng', help ='select the atoms for fengycin',type=str)
parser.add_argument('lipid', help ='select the phosphates for lipid',type=str)
parser.add_argument('first',help = 'start analysis from this iteration', type=int)
parser.add_argument('last',help = 'end analysis at this iteration', type=int)


args = parser.parse_args()

## System specific requirements. Creates Atomic groups for fengycins and lipids.
## 
system = args.dir
source = loos.createSystem(args.model)
feng = loos.selectAtoms(source,args.feng)
phosphates = loos.selectAtoms(source, args.lipid)
increment = 1.0 ## this is the bin width for the progress coordinate histogramming

## This function calculates the centroid of all fengycins and lipids in each frame.
def centroids(feng,lipid, box):
    centroids1 = []
    centroids2 = []
    for i in feng:
        i_cent = i.centroid()
        centroids1.append(i_cent)
    for j in lipid:
        for P in j:
            centroids2.append(P.coords())
    return centroids1, centroids2

## This is the function that histograms the z-height of phosphate w.r.t the distance from fengycins
def dist_hist(cent1, cent2, traj_wt, P_bins, sum_wt):
    for f in cent1:
        for p in cent2:
#            print type(p), type(f)
#            dist = f.distance(p,box)
            dist_vec = f-p
            dist = math.sqrt(dist_vec.x()**2 +dist_vec.y()**2 )
            print dist
            if dist >0.0 and dist < 40.0:
                bin_no = int(dist/increment)
                coord_wt = p.z()* traj_wt
                P_bins[bin_no]+= coord_wt
#                print bin_no
                sum_wt[bin_no]+=traj_wt
#               numpy.append(P_bins,coord_wt,axis=1)
#               numpy.append(sum_wt,traj_wt, axis=1)
#    print P_bins
#    print sum_wt
    return P_bins, sum_wt


root_dir=system + "/"+"traj_segs"
westh5 = system + "/" + "west.h5"
iter = []
outputh5file = system + "/"+ "dist.h5"
iter_sorted = []
### Chosen bins are 40 but you can change it. For this script it is hard wired.
bins_P = [[0.0]]*40
#print "length",len(bins_P[0])
sum_wt = [0.0]*40

## This will sort the iteration directories chronologically.
for x in os.listdir(root_dir):
    if (int(x) >= args.first and int(x)<=args.last):
        iter.append(x)
#print iter
iter_sorted = sorted(iter,key=int)
#print iter_sorted

iterations = numpy.array([])
f=h5py.File(outputh5file,"w")
g= f.create_group("iterations")
## This west.h5 file has the weights and progress coordinates for all trajectories in all iterations. 
west =h5py.File(westh5,"r")
## We have to go through all the output directories. The directories are designed with iteration-bin_no.-trajectory.
## Iteration wise going through each bin for the progress coordinate and then go through each frame of each trajectory inside the bin.
for i in iter_sorted:
    t = int(i)
    iter_dir = root_dir + '/' + i
    seg = []
    for j in os.listdir(iter_dir):
        seg.append(j)
    seg_dir=sorted(seg,key=int)
#    print i
    sg_name = "iterations/iter_00"+ i +"/seg_index"
    sg_data=west[sg_name]
    seg_index = numpy.array(sg_data)
    for c in seg_dir:
#        print i, "\t", c
        l=int(c)
        wt = seg_index[l][0]

        segdcd = root_dir + '/' + i + '/' + c+ '/'+'seg.dcd'
        traj = loos.pyloos.Trajectory(segdcd,source)
        k=0
        for frame in traj:
           feng_mol = feng.splitByMolecule()
           phos = phosphates.splitByMolecule()
           box = source.periodicBox()
           cent1, cent2 =centroids(feng_mol,phos,box)
           new_bins_P, new_sum_wt= dist_hist(cent1,cent2,wt,bins_P, sum_wt)
           bins_P = new_bins_P
           sum_wt = new_sum_wt
           print i, c
#print len(bins_P)
print bins_P
## Histogram the phosphates
avg_phos =[[0.0]]*40
for phos in range(len(bins_P)):
    avg_phos[phos]=bins_P[phos]/sum_wt[phos]

##this can be directly used for plotting on gnuplot
for hist in range(len(avg_phos)):
    for k in range(len(avg_phos[hist])): 
        print hist,"\t",avg_phos[hist][k]


west.close()
