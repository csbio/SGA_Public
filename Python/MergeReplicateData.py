#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# MergeReplicateData.py
#
# Takes an SGA rawdata file as input and determines how many replicates there
#   are for each query. It then spits the data back out renaming replicates to
#   split them up. Explore 2,3,4,5,6. In fact, we might as well explore all 
#   the possibilties within this framework. This scrip
# 
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: May 25, 2012
#
# USAGE:
# MergeReplicateData.py split_param ReplicateList.txt BatchList.txt BigDataFile[.gz]
#
# INPUTS:
#    ReplicateList  head's up so we know which queries to grab
#    BatchList      mapping of query_setid -> batch 
#    BigDataFile contains all of the source data
#
# OUTPUTS: output1 output2
#          Each run splits the replicate set into two complementary groups,
#          then writes one groups (and accompanying batch) to output1 and
#          the other to output2
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

if len(sys.argv) < 4:
    print 'too few arguments (try "-h" for help)'
    exit()

split_param    = sys.argv[1] # ex 1,2,4 -> [1,2,4] [3,5,6,7]
replicate_file = sys.argv[2]
batch_file     = sys.argv[3]
big_data_file  = sys.argv[4]


# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(replicate_file):
    print 'replicate_file "' + replicate_file + '" does not exist'
    exit()
if not os.path.exists(batch_file):
    print 'batch_file "' + batch_file + '" does not exist'
    exit()
if not os.path.exists(big_data_file):
    print 'big_data_file "' + big_data_file + '" does not exist'
    exit()

## Step 1: Read in our our list of queries and the number of replicates to expect
# I'm assumming that replicates are never in the same batch
replicate_queries = {} # query -> num_setids
replicate_fid = open(replicate_file, 'r')
replicate_fid.readline() # header line
for line in replicate_fid:
    line = line.strip().split('\t')
    #if int(line[3]) >= 5:
    if int(line[3]) >= 6:
        if int(line[1]) + int(line[2]) < 2:
            replicate_queries[line[0]] = line[3]
        else:
            sys.stderr.write("problem with query: "+line[0]+"\n")

## Step 2: Hash the batch mapping file
keep_batches = set()
query_setids = {}

batch_fid = open(batch_file, 'r')
for line in batch_fid:
    line = line.strip().split('\t')

    # add this set to the list for this query
    query = '_'.join(line[0].split('_')[:-1])
    qset = line[0].split('_')[-1]

    if query not in query_setids:
        query_setids[query] = set()
    query_setids[query].add(qset)

    if query in replicate_queries:
        keep_batches.add(line[1])


## Step 3: For each query, split the set ids into two groups
group1 = [int(x)-1 for x in split_param.split(',')]
size1 = len(group1)
size2 = {}

setA = {}
for query in replicate_queries:
    setA[query] = [list(query_setids[query])[x] for x in group1]
    size2[query] = len(query_setids[query]) - size1

## Step 4: Iterate through the scorefile
# keep replicate queries, renameing them
# keep anything in keep_batches
# Result can be appended to a short set.

if big_data_file[-3:] == '.gz':
    big_data_fid = fileinput.hook_compressed(big_data_file, 'r')
else:
    big_data_fid = open(big_data_file, 'r')


for line in big_data_fid:
    line = line.strip().split('\t')
    if line[0] in replicate_queries:
        if line[3] in setA[line[0]]:
            line[0] = line[0]+'_A'+str(size1)
        else:
            line[0] = line[0]+'_B'+str(size2[line[0]])
        print('\t'.join(line))
    #elif line[5] in keep_batches:
        #print('\t'.join(line))

