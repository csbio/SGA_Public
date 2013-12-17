#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# RemoveOrphans.py
#
# This script takes an rawdata SGA input, and removes all batches
# with fewer than 3 queries. Result to STDOUT
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: May 10, 2012
#
# USAGE:
# RemoveOrphans.py sgadata > sga.sans_orphan...txt
#
# INPUTS:
#    sgadata: rawdata input file in XX col format
#             [sgadata can be in .gz format]
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

if len(sys.argv) < 2:
    print 'too few arguments (try "-h" for help)'
    exit()

input_file = sys.argv[1]

if input_file[-3:] == '.gz':
    input_fid = fileinput.hook_compressed(input_file, 'r')
else:
    input_fid = open(input_file, 'r')

SGA_QUERY_COL = 0
SGA_BATCH_COL = 5

# Let's make it easy on ourselves and consider only
# queries with one replicate (in a single batch)

# query -> batch
batches = {}
for line in input_fid:
    line = line.split('\t')
    q = line[SGA_QUERY_COL]
    b = line[SGA_BATCH_COL]
    if b not in batches:
        batches[b] = set()
    batches[b].add(q)


keep_batches = set();
for b in batches:
    if len(batches[b]) > 2:
        keep_batches.add(b)


# rewind SGADATA to reprocess
input_fid.seek(0)
for line in input_fid:
    b = line.split('\t')[SGA_BATCH_COL]
    if b in keep_batches:
        print(line),
