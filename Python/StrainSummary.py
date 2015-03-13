#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# StrainSummary.py
#
# This script prints out a summary of number of strain types from an SGA 
#      raw data file (9 col, NOT A SCORE FILE ~12 col!)
# Also: prints the number of unique array plate ids and
#       the number of unique batch ids
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: December 5, 2011
#
# USAGE:
# StrainSummary.py sga_file(.gz)
#
# INPUTS: raw_sga_file or profile_sga_file
#
# OUTPUTS: none (stdout)
#
#############################################################################
"""
    print(HELP_TEXT)
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
    print('too few arguments (try "-h" for help)')
    exit()

sga_file = sys.argv[1]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(sga_file):
    print('sga_file "' + sga_file + '" does not exist')
    exit()

sga_fid = fileinput.hook_compressed(sga_file, 'r')

# default to "raw" input files
QUERY_COL = 0
ARRAY_COL = 1
PLATE_COL = 2
UNIQE_COL = 4
BATCH_COL = 5
if 'raw' in sga_file:
   QUERY_COL = 0
   ARRAY_COL = 1
   PLATE_COL = 2
   UNIQE_COL = 4
   BATCH_COL = 5
elif 'release' in sga_file:
   QUERY_COL = 0
   ARRAY_COL = 2

   PLATE_COL = 0
   UNIQE_COL = 0
   BATCH_COL = 0

# step 1 hash all queries and arrays
query_set = set()
array_set = set()
plate_set = set()
uniqe_hsh = {}
duple_hsh = {}
batch_set = set()
query_strain_count = { 'sn':0, 'tsq':0, 'damp':0, 'dma':0, 'tsa':0, 'y':0, 'unann':0}
array_strain_count = { 'sn':0, 'tsq':0, 'damp':0, 'dma':0, 'tsa':0, 'y':0, 'unann':0}

for line in sga_fid:
   line = line.strip().split()
   query_set.add(line[QUERY_COL])
   array_set.add(line[ARRAY_COL])
   # raw stuff
   if 'raw' in sga_file:
      plate_set.add(line[PLATE_COL])	
      batch_set.add(line[BATCH_COL])
      if line[UNIQE_COL] in uniqe_hsh and uniqe_hsh[line[UNIQE_COL]] != line[QUERY_COL]:
         if line[UNIQE_COL] not in duple_hsh:
            duple_hsh[line[UNIQE_COL]] = set([line[QUERY_COL], uniqe_hsh[line[UNIQE_COL]]])
         else:
            duple_hsh[line[UNIQE_COL]].add(line[QUERY_COL])
      else:
         uniqe_hsh[line[UNIQE_COL]] = line[QUERY_COL]

sga_fid.close()

# step 2 count all possible types
# and store all collab strains for printout
collab = set()
for query in query_set:
    if '_' in query:
        if '_y' in query:
           collab.add(query)
        for strain_type in query_strain_count:
            if strain_type in query:
                query_strain_count[strain_type] = query_strain_count[strain_type]+1

    else:
        query_strain_count['unann'] = query_strain_count['unann']+1
   
for array in array_set:
   if '_' in array:
        for strain_type in array_strain_count:
            if strain_type in array:
                array_strain_count[strain_type] = array_strain_count[strain_type]+1
   else:
        array_strain_count['unann'] = array_strain_count['unann']+1
    
# step 3 print results
type_list = ['sn', 'tsq', 'damp', 'dma', 'tsa', 'y', 'unann']
print('\n' + sga_file + ' summary:\ntype\t',end='')
for strain_type in type_list:
    print(strain_type + '\t',end='')

print('total\nquery\t',end='')
for strain_type in type_list:
    print(str(query_strain_count[strain_type]) + '\t',end='')

print(len(query_set))

print('array\t',end='')
for strain_type in type_list:
    print(str(array_strain_count[strain_type]) + '\t',end='')
print(len(array_set))

if 'raw' in sga_file:
   print('unique array plates: ',end='')
   print(len(plate_set))
   print('unique batchs : ', end='')
   print(len(batch_set))
   print(' ')
   if len(duple_hsh) > 0:
       print('plate collisions:')
       for key in duple_hsh:
           print(key + ':\t' + ', '.join(duple_hsh[key]))

if len(collab) > 0:
   print('Collab Strains:')
   for key in collab:
      print(key) 

