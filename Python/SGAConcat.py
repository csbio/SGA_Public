#!/usr/bin/env python3

# imports and constants
import sys 
import os
import fileinput

def help():
   HELP_TEXT = """ 
#############################################################################
# SGAConcat.py
# 
# This script replicates the contents of SGAfile in the output. Then it 
# concatenates the other data.
#
#   plate numbers are incremented past the max plate number in SGAfile
#   batch ids    are incremented past the max batch number in SGAfile
# 	      Each of these are first assigned unique IDs (from 1)
#        to prevent large gaps in batch ids
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: January 26, 2015
#
# USAGE: 
# SGAConcat.py SGAfile.gz OtherData.gz > outputfile
#  
# INPUTS:
#   SGAfile.gz    gzipped raw sgadata file
#   OtherData.gz  gzipped raw sga data
#
# SGA_file and OtherData may be gzipped 
#
# OUTPUTS: (writes to STDOUT)
#
#############################################################################
"""
   print(HELP_TEXT, file=sys.stderr)
   return

################ MAIN FUNCTION

if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
   help()
   sys.exit()

if len(sys.argv) != 3:
   print('Wrong number of arguments, try "-help"', file=sys.stderr)
   sys.exit()

SGAfile = sys.argv[1]
OtherData = sys.argv[2]

fid_1 = fileinput.hook_compressed(SGAfile, 'r')
fid_2 = fileinput.hook_compressed(OtherData, 'r')

max_plate = 0
max_batch = 0
for line in fid_1:
   if SGAfile[-3:] == '.gz':
      line = line.decode('utf-8').strip()
   else:
      line = line.strip()

   split_line = line.split('\t')
   plate = int(split_line[4])
   batch = int(split_line[5])
   if plate > max_plate:
      max_plate = plate

   if batch > max_batch:
      max_batch = batch

   print(line)

fid_1.close()

print('Other data begins on batch '+ str(max_batch+1) +' plate '+str(max_plate+1), file=sys.stderr)

# rip through file 2 once and map each batch and plate to a new low
other_batches = {}
batches_seen = 0

other_plates = {}
plates_seen = 0

for line in fid_2:

   if OtherData[-3:] == '.gz':
      line = line.decode('utf-8').strip()
   else:
      line = line.strip()

   line = line.split('\t')
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
   line = '\t'.join(line)
   print(line)

fid_2.close()
