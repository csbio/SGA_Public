#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# TransposeClustergram.py
#
# This script removes the dentrogram information from a large cdt file
# so can be opened on a mac.
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: November 15, 2012
#
# USAGE:
# TransposeClustergram.py inputfile > outputfile
#
# INPUTS:
#	inputfile: the cdt to "transpose"
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

all_lines = []
first_col = []
input_fid = open(input_file, 'r')
for line in input_fid:
	line = line.strip().split('\t')
	if 'EWEIGHT' in line:
		line[line.index('EWEIGHT')] = 'GWEIGHT'
	elif 'GWEIGHT' in line:
		line[line.index('GWEIGHT')] = 'EWEIGHT'

	all_lines.append(line)
	
input_fid.close()

# these used to be columns
for row in range(len(all_lines[0])):
	print '\t'.join([line[row] for line in all_lines])
