#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# ReplaceStringsFromMap.py
# See also ReplaceOrfWithCommon.py
#
# This script takes in a mapping file, then changes all instances of 
# strings in the first column of the mapping file to the second. 
# Only complete fields are replaced (expected tabs)
# Good for changing one type of ID to another, or in this case adding
# ID numbers to damp queries with mapping file of the form 
# <Yxxxxx_damp	Yxxxxx_damp1234>
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: February 16, 2012
#
# USAGE:
# ReplaceStringsFromMap.py input_file output_file mapping_file
#
# INPUTS:
#      input_file: file that needs fixing
#      mapping file: cannonical mappings
#
# OUTPUTS:
#      output_file: fixed version of the input_file
#
#############################################################################
"""
    print HELP_TEXT
    return

################ MAIN FUNCTION
# imports and constants
import sys
import os
import re
common_map_delim = '\t'

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) < 4:
    print 'too few arguments (try "-h" for help)'
    exit()

input_file = sys.argv[1]
output_file = sys.argv[2]
common_name_file = sys.argv[3]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(input_file):
    print 'input_file "' + input_file + '" does not exist'
    exit()

if not os.path.exists(common_name_file):
    print 'common_name_file "' + common_name_file + '" does not exist'
    exit()

try:
    output_fid = open(output_file, 'w')
except:
    print 'Error opening output file ' + output_file 
    exit() 

## Step 1: hash common name file
common_name_fid = open(common_name_file, 'r')
common_map = {}
for line in common_name_fid:
    line = line.strip().split(common_map_delim)
    if len(line) > 1:
        common_map[line[0]] = line[1]

common_name_fid.close()

## Step 2: iterate over the input 	
input_fid = open(input_file, 'r')
for line in input_fid:
    line = line.strip().split('\t');
    for i in range(len(line)):
        if line[i] in common_map:
            line[i] = common_map[line[i]]
    line = '\t'.join(line)
    output_fid.write(line)
    output_fid.write('\n');

output_fid.close()
