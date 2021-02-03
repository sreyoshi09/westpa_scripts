#!/usr/bin/env python

import loos
import loos.pyloos
import sys
import numpy
import copy
import os
import h5py


### Data for residue-residue contact map with the color axis showing the how often one residue comes close to another.

####Finds the native contact in the initial structure
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
    for m in range(len(residues1)):
        contacts_res=[]
        for n in range(len(residues2)):
            dist2 = (centroids1[m]-centroids2[n]).length2()
            if dist2 < cutoff2:
                contacts_res.append((m,n))
        contacts.append(contacts_res)
    return contacts


def native_contacts(centroids1, centroids2, cutoff, contacts, box):
####Native contacts for each frame with a smooth cutoff function
    num_contacts = 0.0
    for r1,r2 in contacts:
        dist2 = centroids1[r1].distance2(centroids2[r2], box)
        denom = (dist2/(cutoff**2))**6
        denom = denom + 1
        cont = 1 / denom
        num_contacts += cont
    print  num_contacts/len(contacts)
    return num_contacts/len(contacts)


def centroid_calc(residues1, residues2,box):

    ###Calculates the centroid for all residues
    for res1 in residues1:
        c= res1.centroid()
        centroids1.append(c)

    for res2 in residues2:
        c = res2.centroid()
        centroids2.append(c)
    return centroids1,centroids2

"""#######################################################
Begin main program
####################################################### """


"#System specific information:"

reference_structure_file = sys.argv[1]
model_file = sys.argv[2]
selection_string1 = "resid <=11"
selection_string2 = "resid >= 12 && resid <=88"

cutoff = float(sys.argv[3])

system = str(sys.argv[4])
##Iterations of WEighted ensemble simulation which need to be considered
first = sys.argv[5]
last = sys.argv[6]
### Each trajectory has 11 frames
frame_no = 11
###Using LOOS library create a system from the reference structure or the Xray crystal structure
reference = loos.createSystem(reference_structure_file)
###Select particular residues for which we would like to calculate the residue-residue maps
ref_sel1 = loos.selectAtoms(reference, selection_string1)
### Creates a list with residues from one set of synuclein
ref_sel1_residues = ref_sel1.splitByResidue()
####Select the second set of molecules for which we would like to calculate the residue-residue maps
ref_sel2 = loos.selectAtoms(reference, selection_string2)
### Creates a list with residues from another set of synuclein
ref_sel2_residues = ref_sel2.splitByResidue()


###Call the function which calculates the native contacts for the reference structre
res_natives =  find_native_contacts(ref_sel1_residues,ref_sel2_residues, cutoff)


####Now we do similar things to our system which has a psf file
model =  loos.createSystem(model_file)
sel1 = loos.selectAtoms(model, selection_string1)
sel1_residues = sel1.splitByResidue()
sel2 = loos.selectAtoms(model, selection_string2)
sel2_residues = sel2.splitByResidue()

###File systems for WESTPA
root_dir=system + "/"+"traj_segs"

###This file stores all the weights of each trajectory and the corresponding progress coordinate
westh5 = h5py.File(sys.argv[7],'r')
iter_range = []
iter_sorted = []
###Going throught the file system and selecting the iterations on which we would like to run the analysis
for x in os.listdir(root_dir):
    if (int(x)>=int(first) and int(x)<=int(last)):
        iter_range.append(x)
iter_sorted = sorted(iter_range,key=int)


res_cont = []
bins=numpy.zeros([11,9])
increment = 0.1
i = int(first)
avg_res = numpy.zeros([11,9])
s=0

#####Going into the file system starting from the top most directory
for i in iter_sorted:
        iter_dir = root_dir+ "/" + i
        seg = []
        for j in os.listdir(iter_dir):
                seg.append(j)
        seg_dir = sorted(seg,key=int)
        sg_name = 'iterations/iter_00'+i+"/seg_index"
        ####Get the weights for each trajectory
        sg_data=westh5[sg_name]
        seg_index = numpy.array(sg_data)
        for c in seg_dir:
            ####Go through each trajectory
            seg = iter_dir +"/"+ c + '/seg.dcd'
            traj = loos.pyloos.Trajectory(seg,model)
            k = int(c)
            ####The list that stores the histogram
            res_cont = [None] * 11
            #### Centroids for each set of residues are stored in two lists
            centroids1=[]
            centroids2=[]
            ### Go through each trajectory and native contacts for each frame
            wt = seg_index[k][0]
            for frame in traj:
                box = frame.periodicBox()
                centroids1,centroids2 = centroid_calc(sel1_residues, sel2_residues,box)
                total_natives = reduce(lambda x,y :x+y,res_natives)
                total_cont = native_contacts(centroids1,centroids2 ,cutoff, total_natives,box)
                ### Assign each of these values to the histogram
                for l in range(0,11):
                        print i,c,l
                        res_cont[l] = native_contacts(centroids1,centroids2 ,cutoff, res_natives[l],box)
                        print res_cont[l],res_natives[l]
                        if (total_cont >=0.0 and total_cont <=1.0):
                            bin_no = int((total_cont-0.0)/increment)
                            avg_res[l][bin_no]=avg_res[l][bin_no]+ res_cont[l]#*wt/frame_no
                            bins[l][bin_no]+=1
### Print the values in a text file for gnuplot to plot the data
for i in range(avg_res.shape[0]):
    for j in range(avg_res.shape[1]):
        if bins[i][j]<1:
            avg_res[i][j] = avg_res[i][j]/1.0
        else:
            avg_res[i][j] = avg_res[i][j]/(bins[i][j])
numpy.savetxt('out.txt',avg_res)
