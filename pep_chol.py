#!/usr/bin/env python


###Same as  pep_lip_cont.py but now has the option of incorporating cholesterol in the calculation
import loos
import loos.pyloos
import sys
import numpy
import copy
import os
import h5py
import argparse

parser = argparse.ArgumentParser()

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


root_dir=system + "/"+"traj_segs"
westh5 = system + "/" + "west.h5"
iter = []
outputh5file = system + "/"+ "feng_lipid.h5"
iter_sorted = []
for x in os.listdir(root_dir):    
    if (int(x) >= args.first and int(x)<=args.last):
        iter.append(x)
iter_sorted = sorted(iter,key=int)
iterations = numpy.array([])
f=h5py.File(outputh5file,"w")
g= f.create_group("iterations")
west =h5py.File(westh5,"r")

def lipid_feng(source,feng,lipid,box):
#    feng.reimage()
    individual_feng = feng.splitByMolecule()
    upper_leaflet = loos.selectAtoms(lipid,'resid >= 166 && resid <=352')
    individual_lipid = upper_leaflet.splitByMolecule()
    CAB_fl=0
    
    for i in individual_lipid:
        for j in individual_feng:
            
            cent_1=i.centroid()
            cent_2=j.centroid()
            distance_between1 = cent_1.distance(cent_2,box)
            Sij_fl = 1/(1+(distance_between1/10.0)**6)
            CAB_fl += Sij_fl
            print CAB_fl 
    return CAB_fl




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
           data[l][k]=lipid_feng(source,feng,lipid,box)
           print i,'\t',c,'\t',data[l][k]
           k = k+1
           
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
                                                                
    s=g.create_group(name_ds)
    s.attrs['n_iter']= t
    s.create_dataset("feng-lipid",shape=(len(os.listdir(iter_dir)),6,2),data=final)
west.close()
f.close()
