#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# RawdataReplicateRename.py
#
# This script appends modifies each line in the inputfile by appending the
# batch id to the query field via an '_'. These queries will then be scored
# seperately. 
# 
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: January 27, 2011
#
# USAGE:
# RawdataReplicateRename.py combined_data_file.gz output_file 
#
# INPUTS:
#  combined_data_file.gz
#        raw sga file with all data
#        file can be gzipped (recommended)
#############################################################################
"""
    print HELP_TEXT
    return

################ MAIN FUNCTION
# imports and constants
import sys
import os
import fileinput

SGA_QUERY_COL = 0
SGA_SET_COL   = 3

# Option: override with Batch column instead of SET
# SGA_SET_COL   = 5

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) < 3:
    print 'too few arguments (try "-h" for help)'
    exit()

combined_data_file = sys.argv[1]
output_file = sys.argv[2]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(combined_data_file):
    print 'combined_data_file "' + combined_data_file + '" does not exist'
    exit()
try:
    output_fid = open(output_file, 'w')
except:
    print 'Error opening output file: ' + output_file
    exit() 

## Step 1: Split each line and add _SETID to the first field
if combined_data_file[-3:] == '.gz':
    combined_data_fid = fileinput.hook_compressed(combined_data_file, 'r')
else:
    combined_data_fid = open(combined_data_file, 'r')

SEEN_QUERIES = set()
for line in combined_data_fid:
    line = line.split()
    query = line[SGA_QUERY_COL]
    set   = line[SGA_SET_COL]
    # Add _set to query name
    queryset = query + '_' + set
    line[SGA_QUERY_COL] = queryset
    output_fid.write('\t'.join(line))
    output_fid.write('\n')

combined_data_fid.close()
output_fid.close()
