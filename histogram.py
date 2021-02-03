#!/usr/bin/env python
import loos
import math
import loos.pyloos
import sys
import numpy
import copy
import os
import h5py
from collections import Counter
import argparse


#Take the argument and the quantity of system its taking in
variable = "pcoord"
system = str(sys.argv[1])

root_dir=system + "/"+"traj_segs" 
iter_no = sys.argv[2]
nbins = int(sys.argv[3])

#take data from the h5 file created from distance_old.py script
h5file = system+ sys.argv[4]
#print h5file
f=h5py.File(h5file,"r")

#Calculate the minimum and maximum of the data
min_value = float(sys.argv[5])
max_value = float(sys.argv[6])

#calculate bin size
increment = (max_value-min_value)/nbins
#increment = 0.758585951
bin_mid = numpy.empty([nbins+1])
bin_value=min_value+(increment/2)
for x in range(nbins+1):
    bin_mid[x] = bin_value
    bin_value+=increment



#saving each iteration name from the westpa directory in a list
#iteration = []
#for y in os.listdir(root_dir):
#    iteration.append(y)
#print iteration
#iteration_sorted = sorted(iteration,key=int)

#create histogram for each iteration
#ist_list = []
#histogram = numpy.zeros([len(iteration),nbins])
#print histogram.shape
histogram = numpy.zeros([nbins])
#data_iter = []
#for x in range(len(iteration)):
#    if int(x) == int(iter_no):
dist_name =  "iterations" + "/iter_00000"+ iter_no +"/"+"pcoord"
#rint dist_name
weight_name = "iterations" + "/iter_00000"+ iter_no +"/"+"weights"
#        break

variable = f[dist_name]
data = numpy.array(variable)
#print data
#     data_iter.append(data)
weight = numpy.array(f[weight_name])
#print weight
# take sum of the weights
#        weight_sum = (numpy.sum(weight))*6
#        print weight_sum
#       bin_no = 0
#constructing histogram for one iteration
#        histogram = [0.0]*nbins
#print "i and j"
#print data.shape[0],data.shape[1]
for i in range(data.shape[0]):
#    for j in range(data.shape[1]):
        #print "data",data[i][j][0]
#                for k in range(nbins-1):
#            print "bins", bin_mid[k]
        if (data[i][5][0]>min_value) and (data[i][5][0]<max_value):
            bin_no = int((data[i][5][0]-min_value)/increment)
            #       print "inside if", histogram[k],'\t',weight[i]
            histogram[bin_no]=weight[i][0]+histogram[bin_no]
            print data[i][5][0]#,"\t", bin_no
##storing each histogram of each iteration as a member of a list 
#                hist_list.append(histogram)
# print histogram of a particular iteration
#for l in range(nbins): 
#    print bin_mid[l], "\t",histogram[l]/increment
#E = numpy.empty([histogram.shape[0]])
#for l in range(histogram.shape[0]):
#    if(histogram[l] >0):
#        E[l]= (-0.593)*(math.log(histogram[l]))
#    else:
#        E[l] = float('inf')
#print E
#
#for m in range(nbins):
#    print bin_mid[m],"\t", histogram[m]
