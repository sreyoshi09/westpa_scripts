#!/usr/bin/env python

import loos
import loos.pyloos
import sys
import numpy
import copy
import os
import h5py

west=h5py.File("../dispersed/west.h5","r")
name_subgroup = 'iterations/iter_00000050/seg_index'
subgroup= west[name_subgroup]
seg_index = numpy.array(subgroup)
weights = numpy.zeros( seg_index.shape[0] )
temp = 0
for m in seg_index:
    print m[0]
    

