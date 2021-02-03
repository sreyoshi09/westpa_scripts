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


##Define the numpy arrays to store xyrdfs and weights of the trajectories
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

##Call the xy_rdf tool in loos
def xyrdf(psf,seg,feng,lipid):
    out_file = system + '/xyrdf_check'
    hist_min = str(0)
    hist_max = str(40)
    xyrdf_bins =str(40)
####   Only last frame of the trajectories are stored
    array_call = ['xy_rdf','--skip=5', '--split-mode=by-molecule','--reselect', psf, seg, feng, lipid ,'0','40','40'] 
    subprocess.call(array_call, stdout= open(out_file,'w'))
#   print array_call
    

##Average over the trajectories

feng= args.feng
lipid= args.lipid
rdf_total= []
for p in range(len(traj_hist)):
#### The big numpy array which stores all the xyrdf according the histogramming done on the trajectories

    if len(traj_hist[p])!=0:
        rdf_bin = numpy.zeros([len(traj_hist[p]),40,2])
        for q in range(len(traj_hist[p])):
            xyrdf(psf,traj_hist[p][q],feng,lipid)
            out = system + '/xyrdf_check'
###  Only the xyrdf of the upper leaflet is considered
            rdf_bin[q] = numpy.loadtxt(out, delimiter="\t",usecols=(0,1))

        rdf_red = numpy.zeros([40,2])

        for s in range(rdf_bin.shape[0]):
            for u in range(rdf_bin.shape[1]):
### Averaging the xyrdf by multiplying with corresponding weights over all bins of the xyrdf
                rdf_red[u][1] += rdf_bin[s][u][1]*traj_wt[p][s]/(sum_traj_wt[p])
                rdf_red[u][0] = rdf_bin[0][u][0]

        rdf_total.append(rdf_red) 
    else:
        rdf_red = numpy.zeros([40,2])
        increase = 0.5
        for u in range(0,40):
            rdf_red[u][0] = increase
            increase += 1
        rdf_total.append(rdf_red)

### Print the array so that its a direct input for gnuplot 
increase = 0.1
for d in range(len(rdf_total)):
    for e in range(len(rdf_total[d])):
        print increase,"\t", rdf_total[d][e][0], rdf_total[d][e][1] 
    increase +=0.1

west.close()
