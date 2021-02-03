#!/usr/bin/env python
### This is a sample script to add to w_pdist for getting free energy/probability curves for other quantities apart from the progress coordinates
def getzdistance(n_iter, itergroup):
#### Function is used to read hdf file(dist.h5) containing z-distance between micelle and membrane and then get the probability from w_pdist
    import sys
    import h5py
    import numpy
    distance_file = h5py.File('dist.h5','r')
    iter_zdist = distance_file['iterations/iter_{:08d}/z-distance'.format(n_iter)]
    print iter_zdist
    output_array = numpy.array(iter_zdist)
    return output_array

def get_aggs(n_iter, itergroup):
#### Function is used to read hdf file(agg_1d.h5) containing iterations with just no. of aggregates information in it. w_pdist will output the average probability of different no. of aggregates
    import sys
    import h5py
    import numpy
    aggs_file = h5py.File('agg.h5','r')
    iter_agg = aggs_file['iterations/iter_{:08d}/aggs'.format(n_iter)]
    print iter_agg

    output_array = numpy.array(iter_agg)
    return output_array

def get_fl(n_iter, itergroup):                                          
#### Function is used to read the feng-lipid contacts at each iterations from the feng_lipid.h5 file. This can ultimately give the other h5 file which is plugged into plothist to get nice graphs
    import sys
    import h5py
    import numpy
    fl_file = h5py.File('feng_lipid.h5','r')
    iter_fl = fl_file['iterations/iter_{:08d}/feng-lipid'.format(n_iter)]
    print iter_fl
    output_array = numpy.array(iter_fl)
    return output_array
