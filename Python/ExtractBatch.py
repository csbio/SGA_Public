#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# ExtractBatch.py
#
# This script rips through a rawdata input file and prints out lines with
# the matching batch id or ids
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: December 19, 2011
#
# USAGE:
# ExtractBatch.py SGA_file.gz batch1 batch2 batch3 ... 
#
# INPUTS:
#    SGA_file a rawdata input file (9 col)
#    batch1 ... list of batch ids to extract
#
# OUTPUTS:
#    lines are faithfully reproduced to STDOUT
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

SGA_file = sys.argv[1]
batch_ids = sys.argv[2:]
#BATCH_file = sys.argv[2]
#batch_ids = set()
#b_fid = open(BATCH_file,'r')
#for line in b_fid:
#    batch_ids.add(line.strip())

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(SGA_file):
    print 'SGA_file "' + SGA_file + '" does not exist'
    exit()

if SGA_file[-3:] == '.gz':
    SGA_fid = fileinput.hook_compressed(SGA_file, 'r')
else:
    SGA_fid = open(SGA_file, 'r')


for line in SGA_fid:
    line = line.strip()
    parsed = line.split('\t')
    #if parsed[5] not in batch_ids:
    if parsed[5] in batch_ids:
        print(line)
