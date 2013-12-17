#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# GenerateCoBatchStandard.py
#
# This script takes an rawdata SGA input and generates a list of
# Co-Batch queries to use as a functional standard. If profile
# similarity can predict Co-Batch we're in bad shape.
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: May 10, 2012
#
# USAGE:
# GenerateCoBatchStandard.py sgadata > [ co_batch_std.txt | orphan_list.txt ]
#
# INPUTS:
#    sgadata: rawdata input file in XX col format
#             [sgadata can be in .gz format]
#
# OUTPUTS:
#    list of co-batch pairs is printed to std_out
#    or an orphan list depending on which block is uncommented
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
queries = {}
skip_q = set()
for line in input_fid:
    line = line.split('\t')
    q = line[SGA_QUERY_COL]
    b = line[SGA_BATCH_COL]
    if q not in queries:
        queries[q] = set()
    queries[q].add(b)
    




# now collect the query list for each batch
batches = {}
for q in queries:
    batch_set = queries[q]
    for b in batch_set:
        if b not in batches:
            batches[b] = set()
        batches[b].add(q)

# this prints a list of cobatch pairs
# now print out all pairs for anything with > 1 query
for b in batches:
    Qs = batches[b]
    while len(Qs)>0:
        item1 = Qs.pop()
        for item2 in Qs:
            print(item1+'\t'+item2)


# this prints a list of orphans
#for b in batches:
#    Qs = batches[b]
#    if len(Qs) == 1:
#        print(Qs.pop())
