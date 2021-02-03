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
### also generates weighted average distribution of mers.


parser = argparse.ArgumentParser()
parser.add_argument('dir', help='Path of the system', type=str)
parser.add_argument('model', help='PSF of the system', type=str)
parser.add_argument('feng', help='select atoms for fengycin')
parser.add_argument('first',help = 'start analysis from this iteration', type=int)
parser.add_argument('last',help = 'end analysis at this iteration', type=int)
args = parser.parse_args()


#This function calculates a list of list
# Each comprosing list consists of the fengycin numbers that are part of the list
# for instance if there are 10 fengycins with 3,6,7 forming one aggregate then 1,8,9 form another while 0,2,4,5 form another aggregates
# the output of this function will be ((3,6,7),(1,8,9),(0,2,4,5))
def calculate_aggs_list(molecules_copy,box):
# molecule_copy is a list of each individual fengycin molecules. Each member of the list is an atomic group.
    aggregates = []
####Create an array of pairs of fengycins which form dimers
    pairs = numpy.zeros([15,15])
####Looping over the molecules and calculating the distance between them to determine aggregate sizes 
    for i in range(len(molecules_copy)-1):
        for j in range(i+1,len(molecules_copy)):
            if(molecules_copy[i].contactWith(6.5,molecules_copy[j],box,50)==True) :
                pairs[i][j]=1
    for i in range((pairs.shape[0])-1):
       for j in range(i+1,pairs.shape[1]):
            i_index = []
            j_index =[]
###Depending on two fenguycins in contact we will build the aggregates list
            if pairs[i][j] == 1:
#### If the two fengycin molecules are in contact then what happens
#### if the size of the aggregates are big enough
                for k in range(len(aggregates)):
                    if i in aggregates[k]:
                         i_index.append(k)
                    if j in aggregates[k]:
                         j_index.append(k)
#### to make sure that (1,2) and (2,1) are taken inti account
### if the fengycins in the pairs are not present in aggregates then we make it a new aggregate
                if len(i_index) == 0 and len(j_index)==0:
                     aggregates.append([i,j])
#### If one of the fengycins is already a part of an existing aggregate we append the other member of the pair to this existing aggregate 
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
### Do the same for the other part of the (i,j) pair which is in contact.
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

### If both are parts of aggregates then we merge the two aggregates into one
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
    print "Aggregates are", aggregates
    return aggregates

### From the aggregate list we count the ones whose size is 1 and then the ones whose size is more than 1
def calculate_monomers(list_aggs):
    monomer =[]
    aggs = []
    if len(list_aggs)>0:
        for i in range(len(list_aggs)):
            for j in range(len(list_aggs[i])):
                aggs.append(list_aggs[i][j])
    else:
        for i in range(0,15):
            monomer.append(i)
    if len(list_aggs)>0:
        for i in range(0,15):
            if i not in aggs:
                monomer.append(i)
    print "Monomers are",monomer
    return monomer

### Last function to calculate both monomer numbers and aggregate numbers
def no_of_mers(monomer,aggregates):
    no_aggs = 0.0
    if len(monomer)==0:
       no_aggs = len(aggregates)
    else:
        no_aggs = len(aggregates) + len(monomer)
        print "Aggregates are",no_aggs 
    return no_aggs


### How big are the aggregates. We get a list with how many monomer, dimer, trimers, quartermers... present
def agg_composition(monomer,aggregates):
    mer_comp = numpy.zeros([15])
    for k in range(1,mer_comp.shape[0]):
        for i in range(len(aggregates)):
            if len(aggregates[i])==k+1:
                mer_comp[k]+=k+1
    if len(monomer)>0:
        for i in range(len(monomer)):
            mer_comp[0]+=1
    print "Distribution of fengycins in aggregates",mer_comp
    return mer_comp
    
    



### System specific requirements:
system = args.dir
psffile =args.model
source = loos.createSystem(psffile)
###Selects fengycin molecules only 
fengs = loos.selectAtoms(source,args.feng)
molecules= fengs.splitByMolecule()

#### These are directories associated with WESTPA runs
root_dir=system + "/"+"traj_segs"
#### HDF5 files associated with WESTPA runs
westh5 = system + "/" + "west.h5"
iter = []
agg_1d = system + "/"+ "agg_1d.h5"
agg_2d = system + "/"+ "agg_2d.h5"

####Sort the WESTPA directories according to WESTPA iterations
iter_sorted = []
for x in os.listdir(root_dir):
    if (int(x) >= args.first and int(x)<=args.last):
        iter.append(x)
iter_sorted = sorted(iter,key=int)
iterations = numpy.array([])
###Output HDF5 file with the no. of aggregates and the fengycin-fengycin contacts from the WESTPA run
file_obj_agg_2d=h5py.File(agg_2d,"w")
##Create groups for the 2D file 
group_2d= file_obj_agg_2d.create_group("iterations")
westh5 = system + "/" + "west.h5"
###Input file that we get as output from WESTPA runs
west = h5py.File(westh5,"r")
#####Output file for just the no. of aggregates
file_obj_agg_1d = h5py.File(agg_1d,"w")
###Create group for the 1D HDF5 output file
group_1d= file_obj_agg_1d.create_group("iterations")

###This array will store the composition of the aggregates over trajectory
cluster_comp = []
###Stores the weight of each trajectory for each iteration
traj_wt = []

###Loops through all the iterations calculates no. of aggregates in all framesand all trajectories and stores them in hdf file
for i in iter_sorted:
    iter_dir = root_dir + '/' + i
    t = int(i)
    seg=[]
    name_ds='iter'+'_'+'00'+i
    sg_name = "iterations/iter_00" + i + "/seg_index"
    sg_data = west[sg_name]
    seg_index = numpy.array(sg_data)
### Data file to output no. of aggregates vs probability figure
    data_1d = numpy.zeros([len(os.listdir(iter_dir)),6])
### Data file to output no. of aggregates vs fengycin-fengycin contacts vs probability
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
            print i,"\t", c, "\t", k
            no_aggs = 0.0
            molecules = fengs.splitByMolecule()
            box =  source.periodicBox()
### Calling the different functions
            feng_in_aggs =  calculate_aggs_list(molecules,box)
            feng_in_monomers = calculate_monomers(feng_in_aggs)
            no_of_aggs = no_of_mers(feng_in_monomers,feng_in_aggs)
            agg_comp=  numpy.zeros([15])
            agg_comp = agg_composition(feng_in_monomers,feng_in_aggs) 
            cluster_comp.append(agg_comp)
            traj_wt.append(wt)
            data_1d[l][k] = no_of_aggs
            k = k+1
### How the hdf5 file should look for it to be processed with w_pdist
    name_ds = 'iter' + '_' + '00' + i 
    pcoord_name = 'iterations' + '/'+ name_ds + '/'+ 'pcoord'
    pcoord_data = west[pcoord_name]
    pcoord = numpy.array(pcoord_data)
    final = numpy.empty([pcoord.shape[0],pcoord.shape[1],2])

    for u in range(0,pcoord.shape[0]):
        for v in range(0,pcoord.shape[1]):
            data_2d[u][v][0]=pcoord[u][v][0]
            data_2d[u][v][1]=data_1d[u][v]
    name_subgroup = 'iterations'+'/'+name_ds+'/'+'seg_index'
    s=group_2d.create_group(name_ds)
    s.attrs['n_iter']= t
    s.create_dataset("aggs_2d",shape=(len(os.listdir(iter_dir)),6,2),data=data_2d)
    g= group_1d.create_group(name_ds) 
    g.attrs['n_iter']=t
    g.create_dataset("aggs_1d", shape=(len(os.listdir(iter_dir)),6,1), data=data_1d)

file_obj_agg_2d.close()
file_obj_agg_1d.close()
### Uncomment these lines to get probability of number of aggregates vs probability hdf file
#for i in iter_sorted:
#    for j in 
#    item_dir = root_dir + '/' + i
#    t = int(i)
#    name_ds = 'iter'+ '_'+'00'+i
#    name_subgroup= 'iterations'+'/'+name_ds+'/'+'seg_index'
#    s=group_1d.create_group(name_ds)
#    s.attrs['n_iter']=t
#    s.create_dataset("aggs_1d",shape=(len(os.listdir(iter_dir)),6,1),data=data_1d)


mers_avg = numpy.zeros([15])

sum_traj_wt = sum(traj_wt)
#### If we want to get the mers vs probability data the following part helps in that
for stuff in range(len(cluster_comp)):
    for value in range(len(cluster_comp[stuff])):
        print 'mer',(value+1),'no. of fengs in that mer',cluster_comp[stuff][value],'sum of traj_wt', sum_traj_wt,'traj_wt',traj_wt[stuff], mers_avg[value]
        mers_avg[value] += (cluster_comp[stuff][value]*traj_wt[stuff]*(value+1))/(sum_traj_wt*15)
for mer in range(len(mers_avg)):
    print mers_avg[mer]





#### We will get hdf files which needs to be further processed using w_pdist
## Example of the w_pdist command will be : w_pdist -b '[[-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0],[-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5]]' --construct-dataset add_on.get_aggs --last-iter 440 -W ../../../west.h5 -o na_ff.h5


#####Uncomment this to run the function on single trajectories
#for frame in trajectory:
#    molecules = fengs.splitByMolecule() 
#    box = source.periodicBox()
#    feng_in_aggs =  calculate_aggs_list(molecules,box)
#    print "Aggregates are",#feng_in_aggs
#    feng_in_monomers = calculate_monomers(feng_in_aggs)
#    print  "Monomers are",feng_in_monomers
#    no_of_aggs = no_of_mers(feng_in_monomers,feng_in_aggs)
#    print "No. of aggregates",no_of_aggs
#    agg_comp=  numpy.zeros([15])
#    agg_comp = agg_composition(feng_in_monomers,feng_in_aggs)
#    print "Aggregate composition",agg_comp
