#!/usr/bin/env python

"""
p_val.py

"""

import sys
import numpy
import argparse
from scipy import stats
from math import isnan

# Parse command-line
parser = argparse.ArgumentParser()
parser.add_argument('prefix', help='Prefix output files with this')
parser.add_argument('--files1', help='List of files of condition 1', nargs='+')
parser.add_argument('--files2', help='List of files of condition 2', nargs='+')
args = parser.parse_args()

# Print header
hdr = " " + " ".join(sys.argv)

# Check if there is the same number of files1 and files2
if len(args.files1) != len(args.files2):
    print "There are %d files1 and %d files2." % (len(args.files1), len(args.files2))

# Create arrays to store data
afile = numpy.loadtxt(args.files1[0])
rows = afile.shape[0]
columns = afile.shape[1]

files1 = numpy.zeros((len(args.files1), rows, columns))
files2 = numpy.zeros((len(args.files2), rows, columns))

# Read in data files
for i in range(len(args.files1)):
    print args.files1[i]
    files1[i] = numpy.loadtxt(args.files1[i])
    files2[i] = numpy.loadtxt(args.files2[i])
print "Done reading in data files."

# Compute p-values
p_val = numpy.zeros((rows, columns))

for i in range(rows):
    for j in range(columns):
        t_test = stats.ttest_ind(files1[:,i,j], files2[:,i,j], equal_var=False)
        ap_val = t_test[1]
        if isnan(ap_val):
            ap_val = 1.0
        p_val[i][j] = ap_val

numpy.savetxt(args.prefix + '.p_val', p_val)
print "Done computing p-values."

