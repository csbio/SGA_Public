#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# addStrainIdToOrf.py
#
# This script is for 2 or 3 column with orfs in column 1 (e.g. fitness
# linkage, or bad_array). Replicates t
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: June 30, 2011
#
# USAGE:
# addStrainIdToOrf.py input_file mapping_file output_file
#
# INPUTS:
#           input_file 	data file to fix
#           mapping_file list of ORFs with ids attached 
#
# OUTPUTS:
#           output_file replicate the input file, adds equvilent lines
#                       and repeats lines for orfs with multiple possible ids
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

if len(sys.argv) < 3:
    print 'too few arguments (try "-h" for help)'
    exit()

input_file = sys.argv[1]
mapping_file = sys.argv[2]
output_file = sys.argv[3]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(input_file):
    print 'input_file "' + input_file + '" does not exist'
    exit()

if not os.path.exists(mapping_file):
    print 'mapping_file "' + mapping_file + '" does not exist'
    exit()

try:
    output_fid = open(output_file, 'w')
except:
    print 'Error opening output file ' + output_file
    exit() 


# Step one, hash the mapping file
mapping_fid = open(mapping_file, 'r')
strains = {}
for line in mapping_fid:
    (orf, strain) = line.strip().split('_')
    if orf not in strains:
        strains[orf] = set()

    strains[orf].add(strain)
mapping_fid.close()

# step two, rip through the input file 
input_fid = open(input_file, 'r')

for line in input_fid:
    output_fid.write(line)
    (orf, remainder) = line.split('\t', 1)
    if '_' in orf:
        (orf, junk) = orf.split('_', 1)

    if orf in strains:
        strain_list = strains[orf]
    else:
        continue

    for strain in strain_list:
        output_fid.write(orf + '_' + strain + '\t' + remainder)

output_fid.close()
input_fid.close()
 
