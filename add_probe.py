#! /usr/bin/env python

"""
  add_probe.py

    Add probe(s) chosen randomly from a given library of conformations to a system. 
    The probe is placed within the period box of the system at at least 1.75 Angstroms
    away from any atoms of the "solute" (e.g. protein; protein & detergent micelle).

    Notes: some of the code was adapted from OptimalMembraneGenerator.py

  Usage- 
     python add_probe.py prefix model.pdb solute #probes probes-library

  Example-
     python add_probe.py rhod 1u19_postprocessed.pdb '(segid == "RHOD" || segid == "MEMB") && !hydrogen' 1 lipid_lib/*

"""

import sys
import argparse
import loos
import loos.pyloos
import random
import math

# Parse command line
parser = argparse.ArgumentParser()
parser.add_argument('prefix', help='Prefix output files with this', type=str)
parser.add_argument('model', help ='PDB of the system', type=str)
parser.add_argument('solute', help='LOOS selection for solute', type=str)
parser.add_argument('probes', help='Number of probes to add', type=int)
parser.add_argument('library', help='Path of probe library', type=str)
args = parser.parse_args()

# Header
hdr = "# " + " ".join(sys.argv)

# Create model and select solute
model = loos.createSystem(args.model)
solute = loos.selectAtoms(model, args.solute)
outgrp = model.copy()

# Determine system box size
box = model.boundingBox() # CHARMM-GUI output is not periodic

# Keep track of added probs
added_probes = []

# Probe placement criteria
x_axis = loos.GCoord(1,0,0)
y_axis = loos.GCoord(0,1,0)
z_axis = loos.GCoord(0,0,1)
overlap_dist = 1.75   # the distance that is considered a "contact"
overlap_threshold = 1 # the number of contacts required for a molecule
                      # to be rejected
# Set up lipid library
library = LipidLibrary.LipidLibrary(args.library)
print "%d configurations found in probe library." % library.size()
i=0
# Add probes iteratively
while i < args.probes:

    outside = False
    
    # Pick a probe from the library
    probe = library.pick_structure()
    
    # Put the molecule at the origin
    probe.centerAtOrigin()

    # Perform random rotation about x, y and z axis 
    probe.rotate(x_axis, random.uniform(0.,360.))
    probe.rotate(y_axis, random.uniform(0.,360.))
    probe.rotate(z_axis, random.uniform(0.,360.))

    # Generate new probe coordinates
    x = random.uniform(box[0].x(), box[1].x())
    y = random.uniform(box[0].y(), box[1].y())
    z = random.uniform(box[0].z(), box[1].z())
    
    vec = loos.GCoord(x,y,z)
    
    # Translate probe to new location
    probe.translate(vec)

    # Check that probe is within bounding box
    for atom in probe:
        coords = atom.coords()
        if (coords.x() < box[0].x()) or (coords.x() > box[1].x()) or (coords.y() < box[0].y()) or (coords.y() > box[1].y()) or (coords.z() < box[0].z()) or (coords.z() > box[1].z()):
            outside = True
            break

    if outside: continue

    # Do a bump-check against the solute
    overlap = probe.within(overlap_dist, solute)
    if len(overlap) > overlap_threshold: continue

    # Do a bump-check against previously added probes
    for j in range(len(added_probes)):
        overlap = probe.within(overlap_dist, added_probes[i])
        if len(overlap) > overlap_threshold: break
    if len(overlap) > overlap_threshold: continue

    print "Probe was successfully placed."
    i += 1

    # Remove overlapping waters/ions
    solvent = loos.selectAtoms(outgrp, '!(' + args.solute + ') && !hydrogen')
    clashing_solvent = solvent.within(overlap_dist, probe)
    solvent_residues = solvent.splitByResidue()
    j = 0
    
    for atom in clashing_solvent:
        for residue in solvent_residues:
            if atom in residue:
                outgrp.remove(residue)
                j += 1
                break

    print "%d solvent residues were removed due to clashes." % j
    
    # Add probe to output atomic group
    for atom in probe:
        atom.resname("CHAP")
        atom.segid("FREE")
    outgrp.append(probe)

outgrp.renumber()

# write out a final pdb file
pdbfile = open(args.prefix + ".pdb", "w")
pdb = loos.PDB.fromAtomicGroup(outgrp)
pdbfile.write(str(pdb))
pdbfile.close()
