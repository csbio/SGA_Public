#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# release_summary.py
#
# This script summarizes a release dataset (9 column)
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: October 22, 2012
#
# USAGE:
# profile_summary.py inputfile[.gz] >> output
#
# OUTPUTS:
# filename	queries	sn	damp	tsq	double	other	arrays	dma	tsa	other FG30 FG26 FG
#
#############################################################################
"""
    print HELP_TEXT
    return

def header():
     HEADER_TEXT = ['filename', 'queries', 'sn', 'damp', 'tsq', 'double', 'other', 'arrays', 'dma', 'tsa', 'other', 'Neg', 'Pos', 'Insig', 'NaN', 'FG30', 'FG26', 'TS30', 'TS26', 'SCI']
     print '\t'.join(HEADER_TEXT) 

################ MAIN FUNCTION
# imports and constants
import sys
import os
from fileinput import hook_compressed
from numpy import isnan

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if sys.argv.count('-header') > 0:
    header()
    exit()

if len(sys.argv) < 2:
    print 'too few arguments (try "-h" for help)'
    exit()

input_file= sys.argv[1]
if input_file[-3:] == '.gz':
    input_fid = hook_compressed(input_file, 'r')
else:
    input_fid = open(input_file, 'r')


queries = set()
arrays = set()
counts = [0,0,0,0] # neg pos insig nan
labels = ['FG30', 'FG26', 'TS30', 'TS26', 'SCI']
exp_counts = {}
for i in range(len(labels)):
	exp_counts[labels[i]] = 0

input_fid.readline() # toss header
for line in input_fid:
	line = line.strip().split('\t')
	queries.add(line[0])
	arrays.add(line[2])
	eps = float(line[5])
	pvl = float(line[7])
	exp_counts[line[8]]+=1
	if isnan(eps) or isnan(pvl):
		counts[3]+=1
	elif pvl < 0.05 and eps >=  0.08:
		counts[1]+=1
	elif pvl < 0.05 and eps <= -0.08:
		counts[0]+=1
	else:
		counts[2]+=1


qtypes = ['_sn', '_damp', '_tsq', '+']
qcount = [sum([sub in orf for orf in queries]) for sub in qtypes]
atypes = ['_dma', '_tsa']
acount = [sum([sub in orf for orf in arrays]) for sub in atypes]


# OK, now print out the restults
result = sys.argv[1]+'\t'
result += str(len(queries))+'\t'
result += '\t'.join([str(qcount[i]) for i in range(len(qcount))]) + '\t'
result += str(len(queries) - sum(qcount[:3]))+'\t'
result += str(len(arrays)) +'\t'
result += '\t'.join([str(acount[i]) for i in range(len(acount))])+'\t' 
result += str(len(arrays) - sum(acount))+'\t' 
result += '\t'.join([str(counts[i]) for i in range(len(counts))])+'\t'
result += '\t'.join([str(exp_counts[i]) for i in labels])
print(result)

