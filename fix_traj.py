#!/usr/bin/env python

import os, sys
import os.path
import subprocess
system = "/run/media/sreyoshi/sheriff/feng_1_chol_mod"
root_dir= system+"/"+"traj_segs"
psffile = system +"/" "bstates/feng_1_chol_mod.psf"
output = []
input1 = []
first = int(sys.argv[1])
last = int(sys.argv[2])




for x in os.listdir(root_dir): 
    print int(x)
    if (int(x)>=first and int(x)<=last):
        print x
        iter_dir = root_dir + '/' + x
        for k in os.listdir(iter_dir):
            seg_iter=k
            dcd=iter_dir+'/'+seg_iter+'/seg'
            xtc=iter_dir+'/'+seg_iter+'/seg.xtc'
            output.append(dcd)
            input1.append(xtc)
print output        
for y in range(len(output)):
    output_scalar = output[y]
    output_dcd =  output[y] + ".dcd"
    input_scalar = input1[y]
#    if not os.path.isfile(output_dcd):
    print  output_scalar, psffile, input_scalar
    
    sub  = ['subsetter','-C','segid == "POPC" || segid == "CHOL"', '-P','segid=="F22P"', '--reimage=zealous',output_scalar, psffile, input_scalar]
    subprocess.call(sub)
