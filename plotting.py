import matplotlib.pyplot as plt
import numpy as np

def avg_pdist(hist, midpoints, binbounds):
    plt.xlabel('Fengycin-Fengycin contacts')
    plt.ylabel('Aggregation free energy')

def na_cont(hist, midpoints, binbounds):
    plt.title('Free energy w.r.t no. of aggregates')
    plt.xlabel('Number of Feng-Feng contacts')
    plt.ylabel('No. of aggregates')
    plt.xlim(0,25)
    plt.ylim(0,15)

###This is for the evolution plot
def evol_pdist(hist, midpoints, binbounds):
    plt.title('Evolution Feng-Feng contacts')
    plt.xlabel('No. of Feng-Feng contacts')
    plt.xlim(0,25)

def feng_lipid(hist, midpoints, binbounds):
    plt.title('Free energy w.r.t Feng-Lipid contacts')
    plt.xlabel('feng-feng contacts')
    plt.ylabel('feng-lipid contacts')
    plt.xlim(0,25)
    plt.ylim(0,33)
