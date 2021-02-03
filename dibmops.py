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
parser.add_argument('feng', help ='select the atoms for fengycin',type=str)
parser.add_argument('lipid', help ='select the phosphates for lipid',type=str)
parser.add_argument('first',help = 'start analysis from this iteration', type=int)
parser.add_argument('last',help = 'end analysis at this iteration', type=int)

args = parser.parse_args()


##System details
system = args.dir
psf = args.dir + "/bstates/visual.psf"
root_dir=system + "/"+"traj_segs"
westh5 = system + "/" + "west.h5"

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
traj_hist = [[] for x in xrange(9)]
traj_wt = [[] for x in xrange(9)]
increment = 0.1


##Histogram the trajectories according to the fraction of progress-coordinate
for i in iter_sorted:
    t = int(i)
##  Path of the trajectories
    iter_dir = root_dir + '/' + i
    seg=[]
##  Access the progress coordinates from the WESTPA h5 file and store it in a numpy array
    name_ds='iter'+'_'+'00'+i
    pcoord_name = 'iterations'+'/'+name_ds+'/'+'pcoord'
    pcoords_data = west[pcoord_name]
    pcoord = numpy.array(pcoords_data)
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
        frac_pcoord = pcoord[l][5][0]/25
        if (frac_pcoord> 0.0 and frac_pcoord<=1.0):
            bin_no = int(frac_pcoord/increment)
            seg = root_dir + '/' + i + '/' + c + '/'+'seg.dcd'
### Assigning trajectories to bins according to progress coordinates fraction.
            traj_hist[bin_no].append(seg)
### Assigning weights of the trajectories according to the progress coordinates fraction and stored in the same order in which the trajectories are stored
            traj_wt[bin_no].append(wt)




sum_traj_wt =  [[] for x in xrange(9)]
for l in range(0,9):
    sum_traj_wt[l] = sum(traj_wt[l])


def dibmops(psf,seg,feng,lipid):
     out_file = system + '/dibmop_check'
     array_call = ['dibmops','-R','40','-N','40', feng, lipid, psf, seg]
     subprocess.call(array_call, stdout= open(out_file,'w'))

##Average over the trajectories

feng= args.feng
lipid= args.lipid
#rint feng, lipid
mop_total= []
#print len(traj_hist)
for p in range(len(traj_hist)):
#### If there are any trajectories in the first bin of feng-feng contacts
    if len(traj_hist[p])!=0:
#### The big numpy array which stores all the dibmops according the histogramming done on the trajectories
        mop_bin = numpy.zeros([len(traj_hist[p]),40,2])  
        for q in range(len(traj_hist[p])):
            dibmops(psf,traj_hist[p][q],feng,lipid)
            out = system + '/dibmop_check'
###  Only the dibmop of the upper leaflet is considered
            mop_bin[q] = numpy.loadtxt(out, delimiter="\t",usecols=(0,2))
        mop_red = numpy.zeros([40,2])
        for s in range(mop_bin.shape[0]):
            for u in range(mop_bin.shape[1]):
### Averaging the dibmop by multiplying with corresponding weights over all bins of the dibmop
                mop_red[u][1] += mop_bin[s][u][1]*traj_wt[p][s]/(sum_traj_wt[p])
                mop_red[u][0] = mop_bin[0][u][0]
        mop_total.append(mop_red)
    else:
        mop_red = numpy.zeros([40,2])
        increase = 0.5
        for u in range(0,40):
            mop_red[u][0] = increase
            increase += 1
        mop_total.append(mop_red)

### Print the array so that its a direct input for gnuplot 
increase = 0.1
#print 'len(mop_total)','\t',len(mop_total)
for d in range(len(mop_total)):
#    print 'len(mop_total[d])','\t', len(mop_total[d])
    for e in range(len(mop_total[d])):
        print increase,"\t", mop_total[d][e][0], mop_total[d][e][1] 
    increase +=0.1

west.close()
