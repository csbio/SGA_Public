#!/usr/bin/env python
def help():
    HELP_TEXT = """ 
#############################################################################
# TSAConcat.py
# 
# This script replicates the contents of TSAfile in the output. Then it 
# concatenates the other data.
#
#  Differences from SGAConcat.py
#   plate numbers are NOT CHANGED as they are already unique
#   batch ids     are NOT CHANGED as they already match from t26 to t30
#   set ids      are INCREMENTED to make these screens a "real replicate"
# 	        Each of these are first assigned unique IDs (from 1)
#          to prevent large gaps in batch ids
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: November 10, 2011
#
# USAGE: 
# TSAConcat.py TS26 TS30 outputfile
#  
# INPUTS:
#    outputfile     name of the outputfile to write to
#
# SGA_file and OtherData may be gzipped 
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


if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) != 4:
    print 'Wrong number of arguments, try "-help"'
    exit()


SGAfile = sys.argv[1]
OtherData = sys.argv[2]
outputfile = sys.argv[3]
delim = '\t'

if SGAfile[-3:] == '.gz':
    fid_1 = fileinput.hook_compressed(SGAfile, 'r')
else:
    fid_1 = open(SGAfile, 'r')

if OtherData[-3:] == '.gz':
    fid_2 = fileinput.hook_compressed(OtherData, 'r')
else:
    fid_2 = open(OtherData, 'r')

fid_3 = open(outputfile, 'w')

max_set = 0
for line in fid_1:
    split_line = line.split()
    the_set = int(split_line[3])
    if the_set > max_set:
        max_set = the_set

    fid_3.write(line)

fid_1.close()

# rip through file 2 once and map each set a new low
other_sets = {}
sets_seen = 0

for line in fid_2:
    line = line.split()
    the_set = line[3]
    
    if the_set not in other_sets:
        sets_seen = sets_seen+1
        other_sets[the_set] = sets_seen 

    line[3] = str(max_set + other_sets[the_set])
    line = delim.join(line)
    fid_3.write(line)
    fid_3.write('\n')

fid_2.close()
fid_3.close() 
