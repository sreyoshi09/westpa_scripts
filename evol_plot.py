#!/usr/bin/env python

import loos
import loos.pyloos
import sys
import numpy
import copy
import math
import os
import h5py

pdisth5 =h5py.File(sys.argv[1],'r')
histogram_name = '/histograms'
histogram_data  = pdisth5[histogram_name]
histogram = numpy.array(histogram_data)
#print histogram.shape

first = int(sys.argv[2])
last =  int(sys.argv[3])
subset_array = histogram[first:last:1]
#numpy.seterr(divide='ignore', invalid='ignore')
#log = (-1)*numpy.log(histogram)
#or i in range(subset_array.shape[0]):
#   for j in range(subset_array.shape[1]):
#       for k in range(subset_array.shape[2]):

#neg_log = (-1)*log 
#average = numpy.add.reduce(neg_log,0)/(last-first)
#print average.shape
pcoord =  pdisth5['/midpoints']
n_iter = pdisth5['/n_iter']
#auxiliary = pdisth5['/midpoints_1']

#print pcoord.shape
#print auxiliary.shape
#print pcoord.shape
#print auxiliary.shape

for i in range(histogram.shape[0]):
    for j in range(histogram.shape[1]):
        print i,"\t",pcoord[j],"\t", histogram[i][j]
