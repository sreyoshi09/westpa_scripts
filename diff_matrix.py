#!/usr/bin/env python

import loos
import loos.pyloos
import sys
import numpy
import copy
import os
import argparse
import math

parser = argparse.ArgumentParser()

parser.add_argument('system_1',help='Path of the first file', type=str)
parser.add_argument('system_2',help='Path of the second file', type=str)
args = parser.parse_args()

array_first = numpy.loadtxt(args.system_1,usecols=(0), delimiter = " ")
array_second = numpy.loadtxt(args.system_1,usecols=(1), delimiter = " ")

array_1 = numpy.loadtxt(args.system_1,usecols=(2), delimiter = " ")
array_2 = numpy.loadtxt(args.system_2,usecols=(2), delimiter = " ")

#print array_1

#log_array_1 = numpy.log(array_1)
#log_array_2 = numpy.log(array_2)
#print log_array_1
diff = array_1 - array_2

for i in range(len(diff)):
    print array_first[i], array_second[i], diff[i]
