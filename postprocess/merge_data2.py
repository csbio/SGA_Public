#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# merge_data.py
#
# This script takes in two SGA dataset files and produces a single
# merged output for release. It expects a *.orf to be alongside *.txt scores
#
# WARNING: You have to change a flag manually: see line ~120
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: October 2, 2012
#
# USAGE:
# merge_data.py fg30_datafile fg26_datafile outputfile
# The first dataset overrides the second.  fg30, fg26  or ts26, ts30
#
# INPUT: the input files should be in the "profile_" format
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

if len(sys.argv) < 4:
    print 'too few arguments (try "-h" for help)'
    exit()

score_files = sys.argv[1:3]
out_file  = sys.argv[3]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work

for check in score_files:
   if not os.path.exists(check):
      print 'check "' + check + '" does not exist'
      exit()

orf_files = [i[:-4]+'.orf' for i in score_files]
for check in orf_files:
   if not os.path.exists(check):
      print 'check "' + check + '" does not exist'
      exit()

try:
    out_fid = open(out_file, 'w')
except:
    print 'Error opening output file ' + out_file
    exit() 


# Step one, assemble all the query information.
orf_sets = [set(), set()]
orf_fids = [open(orf_files[i], 'r') for i in range(2)]

# Add any orf from the files that doesn't look like an array
for i in range(2):
   for line in orf_fids[i]:
      if('_dma' not in line and '_tsa' not in line):
         orf_sets[i].add(line.strip())
   orf_fids[i].close()

# FG 30 trumps 26: All 30 these 26
# keep only queries from the first set, ignore them in the second file
fg26_keep = set()
for i in orf_sets[1]:
   if i not in orf_sets[0]:
      fg26_keep.add(i)

score_fids = [open(score_files[i], 'r') for i in range(2)]

# input (from 1) "PROFILE"
# 1 - 4 QoQcAoAc
# 5 - eps
# 6 - std
# 7 - pvl
# 8 - DMF

# TYPE 4
# 1  Query ORF 
# 2  Query gene name
# 3  Array ORF 
# 4  Array gene name
# 5  Double mutant fitness
# 6  Epsilon
# 7  Standard deviation
# 8  P-value
# 9  Experiment
header=[
'Query ORF',
'Query gene name',
'Array ORF',
'Array gene name',
'Double mutant fitness',
'Epsilon',
'Standard deviation',
'P-value',
'Experiment']
out_fid.write('\t'.join(header)+'\n')


#  all of the first (0) and conditionally the second (1)
for line in score_fids[0]:
   line = line.strip().split('\t')
   newline = '\t'.join(line[0:4]) 
   newline = newline + '\t' + line[7] + '\t'
   newline = newline + '\t'.join(line[4:7])
   newline = newline + '\t'+'TS26'+'\n'
   #newline = newline + '\t'+'FG30'+'\n'
   out_fid.write(newline)

for line in score_fids[1]:
   line = line.strip().split('\t')
   if line[0] in fg26_keep:
      newline = '\t'.join(line[0:4]) 
      newline = newline + '\t' + line[7] + '\t'
      newline = newline + '\t'.join(line[4:7])
      newline = newline + '\t'+'TS30'+'\n'
      #newline = newline + '\t'+'FG26'+'\n'
      out_fid.write(newline)

[score_fids[i].close() for i in range(2)]
out_fid.close()
