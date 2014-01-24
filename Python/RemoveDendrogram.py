#!/usr/bin/env python2
def help():
    HELP_TEXT = """
#############################################################################
# RemoveDendrogram.py
#
# This script removes the dentrogram information from a large cdt file
# so can be opened on a mac.
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: July 23, 2012
#
# USAGE:
# RemoveDendrogram.py inputfile > outputfile
#
# INPUTS:
#	inputfile: the cdt to "fix"
#
# OUTPUTS:
#	modified version is written to std_out
#
#############################################################################
"""
    print HELP_TEXT
    return

################ MAIN FUNCTION
# imports and constants
import sys
import os
GID_COL = 0
AID_ROW = 1

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) < 2:
    print 'too few arguments (try "-h" for help)'
    exit()

input_file = sys.argv[1]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(input_file):
    print 'input_file "' + input_file + '" does not exist'
    exit()


input_fid = open(input_file, 'r')
row = 0;

for line in input_fid:
    line = line.strip().split('\t')
    if row < 3 and row > 0:
        line[1] = line[0]
    line = line[0:GID_COL]+line[GID_COL+1:]
    if row != AID_ROW:
        print('\t'.join(line))

    row = row+1
    
