#!/usr/bin/env python

import loos
import loos.pyloos
import sys
import numpy
import copy
import os
import h5py
import argparse
import subprocess

parser = argparse.ArgumentParser()
## This script creates a file with the weighted average molecular order parameter of lipids, distance from the fengycins and the fengycin-fengycin contacts.
##Parse the input from the command line
parser.add_argument('dir',help='Path of the system', type=str)
parser.add_argument('model', help='PSF of the system',type=str)
parser.add_argument('feng', help ='select the atoms for fengycin',type=str)
parser.add_argument('lipid', help ='select the phosphates for lipid',type=str)
parser.add_argument('first',help = 'start analysis from this iteration', type=int)
parser.add_argument('last',help = 'end analysis at this iteration', type=int)

args = parser.parse_args()


##System details
system = args.dir
psf = args.dir + "/bstates/visual.psf"
root_dir=system + "/traj_segs"
westh5 = system + "/west.h5"

##Select the iterations on which averaging is going to be done
iter = []
iter_sorted = []
for x in os.listdir(root_dir):
    if (int(x) >= args.first and int(x)<=args.last):
        iter.append(x)
iter_sorted = sorted(iter,key=int)

### Read the west.h5 file
west =h5py.File(westh5,"r")


##Define the numpy arrays to store dibmops and weights of the trajectories
#traj_hist = [[] for x in xrange(9)]
traj_path =[]
traj_wt = []
increment = 0.1


##Histogram the trajectories according to the fraction of progress-coordinate
for i in iter_sorted:
    t = int(i)
##  Path of the trajectories
    iter_dir = root_dir + '/' + i
    seg=[]
##  Access the progress coordinates from the WESTPA h5 file and store it in a numpy array
    name_ds='iter'+'_'+'00'+i
#    pcoord_name = 'iterations'+'/'+name_ds+'/'+'pcoord'
#    pcoords_data = west[pcoord_name]
#    pcoord = numpy.array(pcoords_data)

## Arrange the iteration directories in chronological order
    for j in os.listdir(iter_dir):
        seg.append(j)
    seg_dir=sorted(seg,key=int)
    sg_name = "iterations/iter_00" + i + "/seg_index"
    sg_data = west[sg_name]
    seg_index = numpy.array(sg_data)
    for c in seg_dir:
        l = int(c)
###Access the weight of the trajectory from the numpy array and store it in the same order as the trajectory list
        wt = seg_index[l][0]
        traj_wt.append(wt)
        traj_name= system +"/traj_segs/" + i +"/"+c + "/seg.dcd"
#        print traj_name
        traj_path.append(traj_name)


###Use this for estimating weighted average
sum_traj_wt = sum(traj_wt)


#for i in range(len(traj_path)):
#    print traj_wt[i],"\t", traj_path[i]

#print sum_traj_wt

feng= args.feng
lipid= args.lipid




### This is the function that calls LOOS to run dibmops on each trajectory which is part of the weighted ensemble run
def dibmops(psf,seg,feng,lipid):
     out_file = system + '/dibmop_check'
     array_call = ['dibmops','-R','40','-N','40', feng, lipid, psf, seg]
     subprocess.call(array_call, stdout= open(out_file,'w'))




mop_avg = numpy.zeros((40,2))
for i in range(len(traj_path)):
    mop = numpy.zeros((40,2))
    dibmops(psf,traj_path[i],feng,lipid)
    out = system +'/dibmop_check'
    mop = numpy.loadtxt(out,delimiter='\t',usecols=(0,2))
#   print rdf
    for j in range(len(mop)):
        mop_avg[j][0]= mop[j][0]
        mop_avg[j][1] += (mop[j][1]*traj_wt[i])/sum_traj_wt



### Final output for being plotted into gnuplot.
for d in range(len(mop_avg)):
        print mop_avg[d][0],"\t",mop_avg[d][1]

west.close()
