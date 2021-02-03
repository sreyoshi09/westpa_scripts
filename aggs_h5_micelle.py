#!/usr/bin/env python

import loos
import loos.pyloos
import sys
import numpy
import copy
import os
import h5py

def no_of_aggs(molecules): 

    molecules_copy = list(molecules)
    clusters=[molecules_copy.pop(0)]
    p = 0
    while(len(molecules_copy)>0):
        non_contacts = []
        flags=0
        for mol in molecules_copy:
            if(mol.contactWith(10,clusters[p],box,10)==True):
                clusters[p].append(mol)
                flag=1
            else:
                non_contacts.append(mol)
        molecules_copy=non_contacts
        if(flags!=1):
            if(len(molecules_copy)>0):
                p=p+1
                clusters.append(molecules_copy.pop(0))
    print  len(clusters)
    return len(clusters)

system = str(sys.argv[1])
psffile =str(sys.argv[2])
source = loos.createSystem(psffile)
fengs = loos.selectAtoms(source,'segid == "F22P" && resname != "TAIL"')
molecules= fengs.splitByMolecule()


root_dir=system + "/"+"traj_segs"
westh5 = system + "/" + "west.h5"
iter = []
outputh5file = system + "/"+ "agg.h5"
iter_sorted = []
for x in os.listdir(root_dir):    
    iter.append(x)
 
iter_sorted = sorted(iter,key=int)
iterations = numpy.array([])
f=h5py.File(outputh5file,"w")
g= f.create_group("iterations")
westh5 = system + "/" + "west.h5"
west = h5py.File(westh5,"r")
#calculat the number of aggregates and put them in a numpy array

for i in iter_sorted:
    print i
    iter_dir = root_dir + '/' + i
    t = int(i)
#     print t
    seg=[]

    data = numpy.empty([len(os.listdir(iter_dir)),6])  
    for j in os.listdir(iter_dir):
        seg.append(j)
    seg_dir=sorted(seg,key=int)
    for c in seg_dir:
        print i, "\t", c
        l=int(c)
        seg = root_dir + '/' + i + '/' + c+ '/'+'seg.dcd'
        print seg
        traj = loos.pyloos.Trajectory(seg,source)
        k=0
        for frame in traj:
            molecules= fengs.splitByMolecule()  
            box = source.periodicBox()
            data[l][k] = no_of_aggs(molecules)
            if(data[l][k] > 15):
                sys.exit("Error message")
            k=k+1
#rieve the pcoord values
    print 'data',data
    print 'data shape',data.shape

    name_ds='iter'+'_'+'00'+i
    pcoord_name = 'iterations'+'/'+name_ds+'/'+'pcoord'
    pcoords_data = west[pcoord_name]
    pcoord = numpy.array(pcoords_data)
    print 'pcoord',pcoord
    print 'pcoord shape',pcoord.shape
    final = numpy.empty([pcoord.shape[0],pcoord.shape[1],2])

    for u in range(0, pcoord.shape[0]):
        for v in range(0,pcoord.shape[1]):
            final[u][v][0]=pcoord[u][v][0]
            final[u][v][1]=data[u][v]

    print 'final',final 
    print 'final',final.shape
    name_subgroup = 'iterations'+'/'+name_ds+'/'+'seg_index'     
    s=g.create_group(name_ds)
    s.attrs['n_iter']= t
    s.create_dataset("aggs",shape=(final.shape),data=final)
f.close()
