#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# AddTSAtoCoords.py
#
# This script takes in an sga file and replicates entries in the coordinate
#  file to satisfy missing annotations
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: November 29, 2011
#
# USAGE:
# AddTSAtoCoords.py sga_file coord_file new_coord_file
#
# INPUTS:
# sga_file: 	must have all the arrays we want to add in column 2
# coord_file:	original file to append
#
# OUTPUTS:
# new_coord_file: same as input with additional entries
#
#############################################################################
"""
    print HELP_TEXT
    return

################ MAIN FUNCTION
# imports and constants
import sys
import os

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) < 4:
    print 'too few arguments (try "-h" for help)'
    exit()

sga_file = sys.argv[1]
coord_file = sys.argv[2]
new_coord_file = sys.argv[3]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(sga_file):
    print 'sga_file "' + sga_file + '" does not exist'
    exit()

if not os.path.exists(coord_file):
    print 'coord_file "' + coord_file + '" does not exist'
    exit()

try:
    new_coord_fid = open(new_coord_file, 'w')
except:
    print 'Error opening output file ' + new_coord_file
    exit() 

# Step one, assemble the set of known things in the coordinate file
# reproduce the input while we're at it
known_coords = {}
coord_fid = open(coord_file, 'r')
for line in coord_fid:
    new_coord_fid.write(line)
    line = line.strip().split('\t', 1)
    known_coords[line[0]] = line[1]

coord_fid.close()

# Step two, find arrays with no coord information
# array is a keyword so i use aarray
new_arrays = set()
sga_fid = open(sga_file, 'r')
for line in sga_fid:
    line = line.strip().split('\t', 2)
    aarray = line[1]
    if aarray not in known_coords:
        new_arrays.add(aarray)

sga_fid.close()

# Step three, print out coord infor for new arrays
for aarray in new_arrays:
    if '_' in aarray:
        (orf, annotation) = aarray.split('_')
        if orf in known_coords:
            new_coord_fid.write(aarray + '\t' + known_coords[orf] + '\n')
        else:
            print('unknowable error' + orf)

    else:
        print('unsplitable error '+ aarray)

new_coord_fid.close()
