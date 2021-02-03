#!/usr/bin/env python
import h5py
import loos
import loos.pyloos
import sys
import numpy
import copy
import os
import argparse


#### This is a script that calculates number of aggregtes of fengycin bound to membrane and stores it in an hdf file
### This will be used by w_pdist in WESTPA to calculate the probability distribution
### also generates weighted average distribution of mers


parser = argparse.ArgumentParser()
parser.add_argument('dir', help='Path of the system', type=str)
parser.add_argument('model', help='PSF of the system', type=str)
parser.add_argument('feng', help='select atoms for fengycin')
parser.add_argument('first',help = 'start analysis from this iteration', type=int)
parser.add_argument('last',help = 'end analysis at this iteration', type=int)
args = parser.parse_args()

### function calculates the no. of aggregates
def no_of_aggs(molecules,box): 
###Makes two lists molecules copy with each fengycin a member of this list and clusters where we keep adding the fengycins which are part of an aggregate
    
    aggregates=[]
    pairs = numpy.zeros([15,15])
    molecules_copy = molecules
    for i in range(len(molecules_copy)-1):
        for j in range(i+1,len(molecules_copy)):
            if(molecules_copy[i].contactWith(6.5,molecules_copy[j],box,50)==True) :
                pairs[i][j]=1
    for i in range(0, pairs.shape[0]-1):
       for j in range(i+1,pairs.shape[1]):
           i_index = []
           j_index =[]
       
           if pairs[i][j] == 1:
               for k in range(len(aggregates)):
                   if i in aggregates[k]:
                        i_index.append(k)
                   if j in aggregates[k]:
                        j_index.append(k)
               if len(i_index) == 0 and len(j_index)==0:
                    aggregates.append([i,j])
               elif len(i_index)>0 and len(j_index)==0:
                    if len(i_index) > 1:
                        newaggs = []
                        newmer = [j]
                        for k in i_index:
                            newmer = list(set(newmer + aggregates[k]))
                        for k in range(len(aggregates)):
                            if k not in i_index:
                                newaggs.append(aggregates[k])
                        newaggs.append(newmer)
                        aggregates = newaggs
                    else:
                        aggregates[i_index[0]].append(j)

               elif (len(j_index)>0 and len(i_index)==0):
                    if len(j_index)>1:
                        newaggs=[]
                        newmer= [i]
                        for k in j_index:
                            newmer = list(set(newmer+aggregates[k]))
                        for k in range(len(aggregates)):
                            if k not in j_index:
                                newaggs.append(aggregates[k])
                        newaggs.append(newmer)
                        aggregates=newaggs
                    else:
                        aggregates[j_index[0]].append(i)
               elif len(i_index)>0 and len(j_index)>0:
                    newaggs = []
                    newmer = []
                    for k in i_index:
                        newmer =list(set(newmer+aggregates[k]))
                    for k in j_index:
                        newmer=list(set(newmer+aggregates[k]))
                    for k in range(len(aggregates)):
                        if ((k not in j_index) and (k not in i_index)):
                            newaggs.append(aggregates[k])
                    newaggs.append(newmer)
                    aggregates= newaggs
    monomer = []
    aggs_list = reduce(lambda x,y : x+y,(aggregates))
    for i in range(0,14):
        if i not in aggs_list:
            monomer.append(i)
    print "Monomers are",monomer
    print "Aggregates are", aggregates

    if len(monomer)==0:
       no_aggs = len(aggregates)
    else:
        no_aggs = len(aggregates) + len(monomer)

    return no_aggs




### System specific requirements:
system = args.dir
psffile =args.model
source = loos.createSystem(psffile)
fengs = loos.selectAtoms(source,args.feng)
molecules= fengs.splitByMolecule()


### Go through the file system again and gain access to west.h5
root_dir=system + "/"+"traj_segs"
westh5 = system + "/" + "west.h5"
iter = []
agg_1d = system + "/"+ "agg_1d.h5"
agg_2d = system + "/"+ "agg_2d.h5"

iter_sorted = []
for x in os.listdir(root_dir):
    if (int(x) >= args.first and int(x)<=args.last):
        iter.append(x)

iter_sorted = sorted(iter,key=int)
iterations = numpy.array([])
file_obj_agg_2d=h5py.File(agg_2d,"w")

group_2d= file_obj_agg_2d.create_group("iterations")
westh5 = system + "/" + "west.h5"
west = h5py.File(westh5,"r")

for i in iter_sorted:

    iter_dir = root_dir + '/' + i
    t = int(i)
    seg=[]
    name_ds='iter'+'_'+'00'+i
    sg_name = "iterations/iter_00" + i + "/seg_index"
    sg_data = west[sg_name]
    seg_index = numpy.array(sg_data)
    data_1d = numpy.zeros([len(os.listdir(iter_dir)),6])
    data_2d = numpy.zeros([len(os.listdir(iter_dir)),6,2])
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
            no_aggs = 0.0
            molecules = fengs.splitByMolecule()
            box = source.periodicBox()
            print i,"\t",c 
            ### Call the no. of aggs function
            no_aggs= no_of_aggs(molecules,box)
            data_1d[l][k] = no_aggs
            #f len(monomer)==0:                                           
            #   no_of_aggs = len(aggregates)
            #else:
            #    no_of_aggs = len(aggregates) + len(monomer)
            #data_1d[l][k]= no_of_aggs
            k = k+1
    name_ds = 'iter'+'_'+'00'+i
    pcoord_name = 'iterations'+'/'+name_ds+'/'+'pcoord'
    pcoords_data =  west[pcoord_name]
    pcoord = numpy.array(pcoords_data)
    final = numpy.empty([pcoord.shape[0],pcoord.shape[1],2])

    for u in range(0, pcoord.shape[0]):
        for v in range(0,pcoord.shape[1]):
            data_2d[u][v][0]=pcoord[u][v][0]
            data_2d[u][v][1]=data_1d[u][v]

    name_subgroup = 'iterations'+'/'+name_ds+'/'+'seg_index'     
    s=group_2d.create_group(name_ds)
    s.attrs['n_iter']= t
    s.create_dataset("aggs_2d",shape=(len(os.listdir(iter_dir)),6,2),data=data_2d)
file_obj_agg_2d.close()

file_obj_agg_1d=h5py.File(agg_1d,"w")
group_1d= file_obj_agg_1d.create_group("iterations")
for i in iter_sorted:
    iter_dir = root_dir + '/' + i
    t = int(i)
    name_ds='iter'+'_'+'00'+i
    name_subgroup = 'iterations'+'/'+name_ds+'/'+'seg_index'
    s=group_1d.create_group(name_ds)
    s.attrs['n_iter']= t
    s.create_dataset("aggs_1d",shape=(len(os.listdir(iter_dir)),6,1),data=data_1d)



file_obj_agg_1d.close()















#####Uncomment this to run the function on single trajectories
#for frame in traj:
#    molecules = fengs.splitByMolecule() 
#    box = source.periodicBox()
#    aggs= no_of_aggs(molecules,box) 
#    print aggs
 
