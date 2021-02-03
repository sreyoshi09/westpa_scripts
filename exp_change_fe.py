#!/usr/bin/env python

import loos
import loos.pyloos
import sys
import numpy
import copy
import os
import h5py
import argparse
###Add functions to calculate native contacts


def find_native_contacts(residues1, residues2, cutoff):
    """ Return a list of residue-residue pairs inside the cutoff distance """

    cutoff2 = cutoff * cutoff
    centroids1 = []
    centroids2 = []

    for res1 in residues1:
        c= res1.centroid()
        centroids1.append(c)

    for res2 in residues2:
        c = res2.centroid()
        centroids2.append(c)

    contacts = []
    for i in range(len(residues1)):
        for j in range(len(residues2)):
            dist2 = (centroids1[i]-centroids2[j]).length2()
            if dist2 < cutoff2:
                contacts.append((i,j))

    print len(contacts)
    return contacts


def native_contacts(residues1, residues2, cutoff2, contacts, box):

    # compute the centroids
    centroids1 = []
    centroids2 = []

    for res1 in residues1:
        c= res1.centroid()
        centroids1.append(c)

    for res2 in residues2:
        c = res2.centroid()
        centroids2.append(c)

    num_contacts = 0.0
    for r1, r2 in contacts:
        dist2 = centroids1[r1].distance2(centroids2[r2], box)
        denom = (dist2 / (cutoff**2))**2
        denom = denom + 1
        cont = 1 / denom
        num_contacts += cont

    return num_contacts/len(contacts)



"""#######################################################
Begin main program
####################################################### """


"#System specific information:"

parser = argparse.ArgumentParser()
parser.add_argument('dir',help='Path of the system', type=str)
parser.add_argument('res_sel_1', help ='select the atoms for fengycin',type=str)
parser.add_argument('res_sel_2', help ='select the phosphates for lipid',type=str)
parser.add_argument('first',help = 'start analysis from this iteration', type=int)
parser.add_argument('last',help = 'end analysis at this iteration', type=int)
parser.add_argument('cutoff', help='cutoff for contacts', type=int)
args = parser.parse_args()









reference_structure_file = args.dir + "/gromacs_config/alpha.gro"
model_file = args.dir + "/bstates/model.psf"
selection_string1 = args.res_sel_1
selection_string2 = args.res_sel_2
cutoff = args.cutoff
system = args.dir

cutoff2 = cutoff*cutoff

reference = loos.createSystem(reference_structure_file)
ref_sel1 = loos.selectAtoms(reference, selection_string1)
ref_sel1_residues = ref_sel1.splitByResidue()
#print ref_sel1_residues 
ref_sel2 = loos.selectAtoms(reference, selection_string2)
ref_sel2_residues = ref_sel2.splitByResidue()
natives = find_native_contacts(ref_sel1_residues,ref_sel2_residues, cutoff)
#print natives
model = loos.createSystem(model_file)
sel1 = loos.selectAtoms(model, selection_string1)
sel1_residues = sel1.splitByResidue()
sel2 = loos.selectAtoms(model, selection_string2)
sel2_residues = sel2.splitByResidue()
outputh5file = system +"/analysis/fe_exp4.h5"

"####WESTPA specific information #####"
root_dir=system + "/"+"traj_segs"
westh5 = system + "/" + "west.h5"
iter = []
iter_sorted = []
for x in os.listdir(root_dir):   
    if (int(x) >= args.first and int(x)<=args.last):
        iter.append(x)
 
iter_sorted = sorted(iter,key=int)
iterations = numpy.array([])
f=h5py.File(outputh5file,"w")
g= f.create_group("iterations")
westh5 = system + "/" + "west.h5"
west = h5py.File(westh5,"r")


"####Calculation over all the segments and iterations####"
for i in iter_sorted:
        iter_dir = root_dir+ '/' + i
        t = int(i)
        seg = []
        data = numpy.empty([len(os.listdir(iter_dir)),2])
        pcoord = numpy.empty([len(os.listdir(iter_dir)),2])
        for j in os.listdir(iter_dir):
                seg.append(j)
        seg_dir = sorted(seg,key=int)
        for c in seg_dir:
            print i, "\t", c
            l = int(c)
            seg = root_dir + '/' + i + '/' + c+ '/'+'seg.dcd'
            print seg
            traj = loos.pyloos.Trajectory(seg,model)
            k=0
            for frame in traj:
                    box = frame.periodicBox()
                    native = native_contacts(sel1_residues, sel2_residues,cutoff2, natives, box)
                    data[l][k]= native
                    k=k+1
        
        name_ds='iter'+'_'+'00'+ i
        pcoord_name = 'iterations' + '/'+ name_ds + '/'+ 'pcoord'
        pcoord_data = west[pcoord_name]
        pcoord = numpy.array(pcoord_data)
        final = numpy.empty([pcoord.shape[0],pcoord.shape[1],2])
        for u in range(0, pcoord.shape[0]):
                for v in range(0,pcoord.shape[1]):
                        final[u][v][0]=pcoord[u][v]
                        final[u][v][1]=data[u][v]

#        print final 
#       print final.shape
        name_subgroup = 'iterations'+'/'+name_ds+'/'+'seg_index'     
        s=g.create_group(name_ds)
        s.attrs['n_iter']= t
        s.create_dataset("res_cont",shape=(len(os.listdir(iter_dir)),2,1),data=data)

f.close()
