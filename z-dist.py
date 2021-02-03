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


parser.add_argument('dir',help='Path of the system', type=str)
parser.add_argument('model', help='PSF of the system',type=str)
parser.add_argument('feng', help ='select the atoms for fengycin',type=str)
parser.add_argument('lipid', help ='select the phosphates for lipid',type=str)
parser.add_argument('first',help = 'start analysis from this iteration', type=int)
parser.add_argument('last',help = 'end analysis at this iteration', type=int)

args = parser.parse_args()

system = args.dir
source = loos.createSystem(args.model)
feng = loos.selectAtoms(source,args.feng)
phosphates = loos.selectAtoms(source, args.lipid)

print system
print args.model
print args.feng




root_dir=system + "/"+"traj_segs"
westh5 = system + "/" + "west.h5"
iter = []
outputh5file = system + "/"+ "dist.h5"
iter_sorted = []
print os.listdir(root_dir)
for x in os.listdir(root_dir):
    if (int(x) >= args.first and int(x)<=args.last):
        iter.append(x)
print iter
iter_sorted = sorted(iter,key=int)
print iter_sorted
iterations = numpy.array([])
f=h5py.File(outputh5file,"w")
g= f.create_group("iterations")
west =h5py.File(westh5,"r")

def z_distance(source,feng,phosphates,box):
    lower_leaflet = loos.AtomicGroup()
    upper_leaflet = loos.AtomicGroup()
    for mol in phosphates:
        if (mol.coords()).z() <= 0.0:
            lower_leaflet.append(mol)
        if (mol.coords()).z() >0.0:
            upper_leaflet.append(mol)

    feng_centroid=feng.centroid()
    if feng_centroid.z() <= 0.0:
        lipid_centroid = lower_leaflet.centroid()
        diff = abs(lipid_centroid.z()-feng_centroid.z())
        print "bottom","\t",diff
    if feng_centroid.z()>0.0:

        lipid_centroid = upper_leaflet.centroid()
        diff = abs(lipid_centroid.z()-feng_centroid.z())
        print "top","\t",diff
    return diff

for i in iter_sorted:
    t = int(i)
    iter_dir = root_dir + '/' + i
    seg=[]

    data = numpy.empty([len(os.listdir(iter_dir)),7])
#    name_ds='iter'+'_'+i
    for j in os.listdir(iter_dir):
        seg.append(j)
    seg_dir=sorted(seg,key=int)
    for c in seg_dir:
        print i, "\t", c
        l=int(c)
        seg = root_dir + '/' + i + '/' + c+ '/'+'seg.dcd'
        traj = loos.pyloos.Trajectory(seg,source)
        k=0
        for frame in traj:
           box = source.periodicBox()
           data[l][k]=z_distance(source,feng,phosphates,box)
           k = k+1
    name_ds='iter'+'_'+'00'+i
    name_subgroup = 'iterations'+'/'+name_ds+'/'+'seg_index'

    s=g.create_group(name_ds)
    s.attrs['n_iter']= t
    s.create_dataset("z-distance",shape=(len(os.listdir(iter_dir)),6,1),data=data)
#    seg_index

#    print seg_index
west.close()
f.close()
