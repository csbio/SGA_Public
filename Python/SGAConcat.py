#!/usr/bin/env python
def help():
    HELP_TEXT = """ 
#############################################################################
# SGAConcat.py
# 
# This script replicates the contents of SGAfile in the output. Then it 
# concatenates the other data.
#
#   plate numbers are incremented past the max plate number in SGAfile
#   batch ids     are incremented past the max batch number in SGAfile
# 	        Each of these are first assigned unique IDs (from 1)
#          to prevent large gaps in batch ids
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: November 10, 2011
#
# USAGE: 
# SGAConcat.py SGAfile.gz OtherData.gz outputfile
#  
# INPUTS:
#    SGAfile.gz     gzipped raw sgadata file
#    OtherData.gz  gzipped raw sga data
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

    fid_3.write(line)

fid_1.close()

print('Other data begins on batch '+ str(max_batch+1) +' plate '+str(max_plate+1))

# rip through file 2 once and map each batch and plate to a new low
other_batches = {}
batches_seen = 0

other_plates = {}
plates_seen = 0

for line in fid_2:
    line = line.split()
    plate = line[4]
    batch = line[5]
    
    # assign unique ids to other_file batches and plates on the fly
    if plate not in other_plates:
        plates_seen = plates_seen+1
        other_plates[plate] = plates_seen
    if batch not in other_batches:
        batches_seen = batches_seen+1
        other_batches[batch] = batches_seen

    line[4] = str(max_plate +  other_plates[plate])
    line[5] = str(max_batch + other_batches[batch])
    line = delim.join(line)
    fid_3.write(line)
    fid_3.write('\n')

fid_2.close()
fid_3.close() 
