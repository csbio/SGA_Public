#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# CountCompleteReplicates.py
# This script reports how many "complete" replicates of each SGA query exist
# in the input dataset.
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: August 23, 2011
#
# USAGE: [-h | -help | --help]
# CountCompleteReplicates.py array_size input_file[.gz] > output_file
#
# INPUTS:
# array_size - number of plate to count as "complete"
# input_file - SGA raw input
#            - may be gzipped
#
# OUTPUTS:
#  STDOUT - list [query_orf, small_sets partial_sets complete_sets]]
#  
#
#############################################################################
"""
    print(HELP_TEXT)
    return

################ MAIN FUNCTION
# imports and constants
import sys
import os
import fileinput
from math import floor

SGA_QUERY_COL = 0
SGA_ARRAY_COL = 2 # array plate number
SGA_SETID_COL = 3 # set


## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) < 3:
    print('too few arguments (try "-h" for help)')
    exit()

ARRAY_PLATE_SET_SIZE = int(sys.argv[1])
ARRAY_PLATE_SET_PART = int(floor(ARRAY_PLATE_SET_SIZE / 2))
input_file = sys.argv[2]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(input_file):
    print('input_file"' + input_file + '" does not exist')
    exit()

if input_file[-3:] == '.gz':
    input_fid = fileinput.hook_compressed(input_file, 'r')
else:
    input_fid = open(input_file, 'r')


# Step 1: hash everything...
queries = {}
set_data = {}
for line in input_fid:
    line = line.split()
    query = line[SGA_QUERY_COL]
    setid = line[SGA_SETID_COL]
    array = line[SGA_ARRAY_COL]

    # for each query, keep track of this setid
    if query in queries:
        queries[query].add(setid)
    else:
        queries[query] = set()
        queries[query].add(setid)

    # for each unique combination of query and set we want to keep track of all array plates
    if query+setid in set_data:
        set_data[query+setid].add(array)
    else:
        set_data[query+setid] = set()
        set_data[query+setid].add(array)
    
# Step 2: for every query iterate over every set. each set needs xxx to be complete

print('Orf' + '\t' + 'partial (>0)' + '\t' + 'partial (>' + str(ARRAY_PLATE_SET_PART) + ')' + '\t' + 'full ('+ str(ARRAY_PLATE_SET_SIZE) + ')')
for query in queries:
    complete = 0
    partial = 0
    incomplete = 0
    for setid in queries[query]:
        plates = len(set_data[query+setid])
        if plates > ARRAY_PLATE_SET_SIZE:
            print('error: too many plates '+ query + ' ' + setid)

        elif plates == ARRAY_PLATE_SET_SIZE:
            complete += 1

        elif plates > ARRAY_PLATE_SET_PART:
            partial += 1

        else:
            incomplete += 1

    print(query + '\t' + str(incomplete) + '\t' + str(partial) + '\t' + str(complete))

      

