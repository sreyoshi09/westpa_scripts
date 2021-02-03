#!/usr/bin/env python
import subprocess
import sys

file=open('istate_map.txt','r+').readlines()

for line in iter(file):
    line_array = line.split(" ")
#    print line_array
    copy_list = ["cp","-r"]
#   segs = int(seg)
    print line_array[0], line_array[1], line_array[2]
    segs= int(line_array[1]) 
    copy_string = "../dispersed/traj_segs/0000"+ line_array[0] +"/{number:06}".format(number=segs)+"/seg.gro"
    paste_string = "bstates/"+ line_array[2].rstrip()+ ".gro"
    copy_list = ["cp","-r", copy_string, paste_string]
    print copy_list
    subprocess.call(copy_list)
