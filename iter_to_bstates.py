#!/usr/bin/env python

import loos
import loos.pyloos
import sys
import numpy
import copy
import os
import h5py
import subprocess



##########Code to copy gro files from segments in traj_seg directory of WESTPA into bstates 
system = str(sys.argv[1])

root_dir = system + "/"+"traj_segs/000"+sys.argv[2]

for x in os.listdir(root_dir):
    seg = root_dir+'/'+x+'/seg.gro'
    command = ['cp',seg,"8mer_new/bstates/"+x+".gro"]
    subprocess.call(command)
