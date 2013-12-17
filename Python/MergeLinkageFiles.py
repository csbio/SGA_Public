#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# MergeLinkageFiles.py
#
# This script merges two linkage files. They may come from different arrays
# (at different chromosomal densities) and we take a conservative approach
# by combining the windows [min(left ends) :  max(right ends)]
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: April 12, 2013
#
# USAGE:
# MergeLinkageFiles.py l_file1 l_file2 > result
#
# INPUTS:
# 	two linkage files to merge
#
# OUTPUTS:
#	Merged file written to STDOUT
#
#############################################################################
"""
    print HELP_TEXT
    return

################ MAIN FUNCTION
# imports and constants
import sys

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) < 3:
    print 'too few arguments (try "-h" for help)'
    exit()

l_file1 = sys.argv[1]
l_file2 = sys.argv[2]

## Step 1: Hash the input data1
linkage_data = {}

l_fid1 = open(l_file1, 'r')
for line in l_fid1:
	line = line.strip().split('\t')
	query = line[0].split('_')[0]
	if query not in linkage_data:
		linkage_data[query] = line[1:]
	else:
		linkage_data[query][0] = str(min(int(linkage_data[query][0]), int(line[1])))
		linkage_data[query][1] = str(max(int(linkage_data[query][1]), int(line[2])))

l_fid1.close()

## Do the same for the second file
l_fid2 = open(l_file2, 'r')
for line in l_fid2:
	line = line.strip().split('\t')
	query = line[0].split('_')[0]
	if query not in linkage_data:
		linkage_data[query] = line[1:]
	else:
		linkage_data[query][0] = str(min(int(linkage_data[query][0]), int(line[1])))
		linkage_data[query][1] = str(max(int(linkage_data[query][1]), int(line[2])))

l_fid2.close()

## Step 2: print out the results
for key in linkage_data:
	print('\t'.join([key, linkage_data[key][0], linkage_data[key][1]]))

