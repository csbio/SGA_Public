#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# ReplaceOrfWithCommon.py
#
# This script takes in a file to change and a common name reference.
# It replaces all instances of orfs (or orfs_junk) with respective common 
# names (or common_junk).
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: November 07, 2011
#
# USAGE:
# ReplaceOrfWithCommon.py input_file output_file common_name_file
#
# Here's a common_name_file:
# /project/csbio/benjamin/Data/Master_Common_Ref_SGD.txt
#
# INPUTS:
#      input_file: file that needs fixing
#      common_name_file: cannonical mappings
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
orf_pattern = re.compile('Y[A-P][L|R][0-9]{3,3}[W|C](?:\-[A-Z]+)?')
input_fid = open(input_file, 'r')
for line in input_fid:
    orfs = re.findall(orf_pattern, line)
    if orfs != []:
        for orf in orfs: 
            if orf in common_map:
                line = line.replace(orf, common_map[orf])

    output_fid.write(line)

output_fid.close()
