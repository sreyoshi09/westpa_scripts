#!/usr/bin/env python

import sys
import numpy
import os

root_dir = str(sys.argv[1])
cum = 0
for i in os.listdir(root_dir):
    iter_dir =  root_dir + '/' + i
    for j in os.listdir(iter_dir):
        cum += 5 
print cum
