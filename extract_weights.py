#!/usr/bin/python

import h5py, numpy, sys
import subprocess
infile = numpy.loadtxt(sys.argv[1], usecols = (0, 2))
itr = sys.argv[2]
sum_wt =0
for iteration, weights in infile[1:]:
    sum_wt += weights
sum_frames =  itr*6
new_frame = 0
for iteration, weights in infile[1:]:
    for frame in range(0,6):
        new_frame = new_frame +1
        new_wt = weights/sum_wt
        print int(new_frame), new_wt

