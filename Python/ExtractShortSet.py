#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# ExtractShortSet.py
#
# This script assembles a small set of various sga queries. We need a small
# set to reduce the iteration time of data scoring. Double and Double data
# must be combined first into one file. First extracts all queries with
# a '+', then adds all SM queries seen in the double mutants as well as any
# other queries in that batch. 
# 
# DEACTIVATED: Then grabs another N_rand random BATCH ids, 
# and add all queries from those batches.
#
# WT data is included, though not explicitly. "undefined+YDL227C" is picked 
# up as a double mutant query and "undefined" is then selected as a single
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: May 12, 2011
# Tested on: Python 2.5.2
#
# USAGE:
# ExtractShortSet.py combined_data_file.gz output_file
#
# INPUTS:
#  combined_data_file.gz
#        raw sga file with all data
#        file can be gzipped (recommended)
#
# You can set SAVE_DEBUG internally to dump pickled data structures
# to a binary output file for analysis and debugging
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
import random
import cPickle

SGA_QUERY_COL = 0
SGA_BATCH_COL = 5
N_rand = 5
SAVE_DEBUG = False

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) < 2:
    print 'too few arguments (try "-h" for help)'
    exit()

combined_data_file = sys.argv[1]
output_file = sys.argv[2]

# Now ensure that these all exist and we're allowed to write the output
if not os.path.exists(combined_data_file):
    print 'combined_data_file "' + combined_data_file + '" does not exist'
    exit()
try:
    output_fid = open(output_file, 'w')
except:
    print 'Error opening output file: ' + output_file
    exit() 


## Step 1: Rip through the file once an get a comprehensive list of batches and queries
# Store a DICT of queries with a SET of batches as values
# Also store a complete set of BATCHES for quick reference
BatchToQuery = {}
QueryToBatch = {}

DoubleQueries = set()
OtherBatches = set()

PrintQueries = set()
PrintBatches = set()

if combined_data_file[-3:] == '.gz':
    combined_data_fid = fileinput.hook_compressed(combined_data_file, 'r')
else:
    combined_data_fid = open(combined_data_file, 'r')

for line in combined_data_fid:
    line = line.split()
    query = line[SGA_QUERY_COL]
    batch = line[SGA_BATCH_COL]
    if query not in QueryToBatch:
        QueryToBatch[query] = set()

    if batch not in BatchToQuery:
        BatchToQuery[batch] = set()

    QueryToBatch[query].add(batch)
    BatchToQuery[batch].add(query)

    if query.find('+') > 0:
        DoubleQueries.add(query)
        PrintQueries.add(query)
        PrintBatches.add(batch)
    else:
        OtherBatches.add(batch)
combined_data_fid.close()


### Step 2: Go DoubleQueries and learn all the single queries we need
# and their batches
for key in DoubleQueries:
    (q1, q2) = key.split('+')
    if q1 in QueryToBatch:
        PrintQueries.add(q1)
        for batch in QueryToBatch[q1]:
            PrintBatches.add(batch)
            OtherBatches.discard(batch)

    if q2 in QueryToBatch:
        PrintQueries.add(q2)
        for batch in QueryToBatch[q2]:
            PrintBatches.add(batch)
            OtherBatches.discard(batch)

## Step 3: Pick N random batches from remaining OTHER batches
#AB_list = list(OtherBatches)
#random.shuffle(AB_list)
#RandomBatches = set(AB_list[:N_rand])
#PrintBatches = PrintBatches.union(RandomBatches)

## Step 4: Get All the Queries in PrintBatches for complete screen printing
for batch in PrintBatches:
    for query in BatchToQuery[batch]:
        PrintQueries.add(query)


## Step 5: Go through the file AGAIN(2) and print anything with a batch in the print set

if combined_data_file[-3:] == '.gz':
    combined_data_fid = fileinput.hook_compressed(combined_data_file, 'r')
else:
    combined_data_fid = open(combined_data_file, 'r')

seen_lines = 0
printed_lines = 0

for line in combined_data_fid:
    seen_lines = seen_lines + 1
    line = line.split()
    query = line[SGA_QUERY_COL]
    batch = line[SGA_BATCH_COL]
    if (batch in PrintBatches) or (query in PrintQueries):
        printed_lines = printed_lines + 1
        output_fid.write('\t'.join(line))
        output_fid.write('\n');

output_fid.close()

## Print a summary
print 'TotalLn: ' + str(seen_lines) 
print 'Printed: ' + str(printed_lines) 
print 'Queries: ' + str(len(PrintQueries))
print 'Batches: ' + str(len(PrintBatches))

save_stuff = [PrintBatches, OtherBatches, PrintQueries, QueryToBatch, BatchToQuery, RandomBatches, DoubleQueries]
comment = """

Load instructions
fid = open(picke_file, 'rb')
save_stuff = cPickle.load(fid)
fid.close()
(PrintBatches, OtherBatches, PrintQueries, QueryToBatch, BatchToQuery, RandomBatches, DoubleQueries) = save_stuff


"""
if SAVE_DEBUG:
    save_fid = open(output_file + '.save.bin', 'wb')
    cPickle.dump(save_stuff, save_fid)
    save_fid.close()
