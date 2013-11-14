#!/usr/bin/env python


def help():

    HELP_TEXT = """
#############################################################################
# merge_arrays.py
#
# This script takes in two merged SGA datasets and produces a single
# merged output for release. Only the first occurrance of any Q/A
# combo is printed.
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: September 27, 2013
#
# USAGE:
# merge_arrays.py fg_merged ts_merged > outputfile
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

score_files = sys.argv[1:3]
score_fids = ['' for i in range(2)]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work

for i in range(2):
    if not os.path.exists(score_files[i]):
        print 'check "' + score_files[i] + '" does not exist'
        exit()

for i in range(2):
    score_fids[i] = fileinput.hook_compressed(score_files[i], 'r')

# Print any combination we haven't seen, as it appears
seen_interactions = {}

for line in score_fids[0]:
    parsed = line.strip().split('\t')
    query = parsed[0]
    array = parsed[2]
    if query not in seen_interactions:
        seen_interactions[query] = set()
    if array not in seen_interactions[query]:
        seen_interactions[query].add(array)
        print(line),

for line in score_fids[1]:
    parsed = line.strip().split('\t')
    query = parsed[0]
    array = parsed[2]
    if query not in seen_interactions:
        seen_interactions[query] = set()
    if array not in seen_interactions[query]:
        seen_interactions[query].add(array)
        print(line),

[score_fids[i].close() for i in range(2)]
