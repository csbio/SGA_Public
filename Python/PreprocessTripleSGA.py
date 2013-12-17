#!/usr/bin/env python
def help():
    HELP_TEXT = """ 
#############################################################################
# PreprocessTripleSGA.py
# 
# This script replicates the contents of SGAfile in the output. Then it 
# concatenates the triple data, making several changes along the way:

#   plate numbers are incremented past the max plate number in SGAfile
#   batch ids     are incremented past the max batch number in SGAfile
#             Each are first decremented by the min, to prevent big
#             gaps from data processed twice
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: June 14, 2012
#
# USAGE: 
# PreprocessTripleSGA.py SGAfile.gz TripleData.gz > outputfile
#  
# INPUTS:
#    SGAfile.gz     gzipped raw sgadata file
#    TripleData.gz  gzipped raw Triple sga data
#    outputfile     script now writes to stdout
#
# SGA_file and TripleData may be gzipped (recommended) or not
#
#############################################################################
"""
    sys.stderr.write(HELP_TEXT)
    return

################ MAIN FUNCTION
# imports and constants
import sys 
import os
import fileinput


if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) != 3:
    sys.stderr.write('Wrong number of arguments, try "-help"')
    exit()


SGAfile = sys.argv[1]
TripleData = sys.argv[2]
delim = '\t'

if SGAfile[-3:] == '.gz':
    fid_1 = fileinput.hook_compressed(SGAfile, 'r')
else:
    fid_1 = open(SGAfile, 'r')

if TripleData[-3:] == '.gz':
    fid_2 = fileinput.hook_compressed(TripleData, 'r')
else:
    fid_2 = open(TripleData, 'r')

max_plate = 0
max_batch = 0
for line in fid_1:
    split_line = line.split()
    plate = int(split_line[4])
    batch = int(split_line[5])
    if plate > max_plate:
        max_plate = plate

    if batch > max_batch:
        max_batch = batch

    print(line),

fid_1.close()

sys.stderr.write('Triple data begins on batch '+ str(max_batch+1) +' plate '+str(max_plate+1))

# Changes to make in file 2:
# increment unique batch ids and plate ids
min_plate = sys.maxint
min_batch = sys.maxint
for line in fid_2:
    line = line.split()
    min_plate = min(min_plate, int(line[4]))
    min_batch = min(min_batch, int(line[5]))

fid_2.seek(0)


for line in fid_2:
    line = line.split()
    line[4] = str(int(line[4]) + max_plate)
    line[5] = str(int(line[5]) + max_batch)
    line = delim.join(line)
    print(line)

fid_2.close()
