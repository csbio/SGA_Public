#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# PurgeDuplicatePlates.py
#
# This script enforces our assumption that unique plate ids are unique.
# It holds the ENTIRE input in memory, if you want to run it on a large
# dataset, re-implement in 2-pass mode, or an on-the-fly remapping of
# more than just the affected plates.
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: March 22, 2012
#
# USAGE:
# PurgeDuplicatePlates.py inputfile[.gz] outputfile
#
#############################################################################
"""
    print HELP_TEXT
    return

################ MAIN FUNCTION
# imports and constants
import sys
import os
import fileinput

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) < 3:
    print 'too few arguments (try "-h" for help)'
    exit()

input_file  = sys.argv[1]
output_file = sys.argv[2]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(input_file):
    print 'input_file "' + input_file + '" does not exist'
    exit()

try:
    output_fid = open(output_file, 'w')
except:
    print 'Error opening output file ' + output_file
    exit() 

if input_file[-3:] == '.gz':
    input_fid = fileinput.hook_compressed(input_file, 'r')
else:
    input_fid = open(input_file, 'r')


IDtoQUERY = {}
input_array = []
plates_to_remap = set();
max_plate_id = 0;

for line in input_fid:
    line = line.split()
    input_array.append(line)
    query = line[0]
    plate = line[4]
    if int(plate) > max_plate_id:
        max_plate_id = int(plate)

    if plate in IDtoQUERY and IDtoQUERY[plate] == query:
        continue
    else:
        plates_to_remap.add(query+plate)
 
new_ids = {}
for line in input_array:
    query = line[0]
    plate = line[4]
    if query+plate in plates_to_remap:
        if query+plate in new_ids:
            line[4] = str(new_ids[query+plate])
        else:
            max_plate_id += 1;
            new_ids[query+plate] = max_plate_id;
            line[4] = str(new_ids[query+plate])
    output_fid.write('\t'.join(line)+'\n')
