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

##Parse the input from the command line
parser.add_argument('dir',help='Path of the system', type=str)
parser.add_argument('model', help='PSF of the system',type=str)
#parser.add_argument('feng', help ='select the atoms for fengycin',type=str)
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


##Define the numpy arrays to store xyrdfs and weights of the trajectories
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



sum_traj_wt = sum(traj_wt)
    

#for i in range(len(traj_path)):
#    print traj_wt[i],"\t", traj_path[i]

#print sum_traj_wt

#feng= args.feng
lipid= args.lipid

def densitydist(psf,seg,lipid):   #feng,lipid):
    out_file = system + '/density_check'
    hist_min = str(-50)
    hist_max = str(50)
    dens_bins =str(400)
    array_call = ['density-dist','--type','mass','--',psf,seg, '-50','50','400',lipid]    #feng,lipid]
    subprocess.call(array_call, stdout= open(out_file,'w'))

dens_avg = numpy.zeros((400,2))
for i in range(len(traj_path)):
    dens = numpy.zeros((400,2))
    densitydist(psf,traj_path[i],lipid)   #feng,lipid)
    out = system + '/density_check'
    dens=  numpy.loadtxt(out,delimiter='\t',usecols=(0,2))
    print dens
    for j in range(len(dens)):
        dens_avg[j][0]= dens[j][0]
        dens_avg[j][1] += (dens[j][1]*traj_wt[i])/sum_traj_wt



for d in range(len(dens_avg)):
        print dens_avg[d][0],"\t",dens_avg[d][1]

west.close()
