#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# AddT30EntriesToFitnessLinkage.py                                           
#                                                                            
# As of today (110620) the linkage and fitness files do not have any info
# for the T30 alleles. This script will duplicate entries with appended
# names. This approach leaves a place to put T30 specific info if they turn
# out to have different linkage regions, or single mutant fitnesses.
#                                                                            
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)                          
# Revision: June 20, 2011
# Tested on: Python 2.6.5                                                    
#                                                                            
# USAGE:                                                                     
# AddT30EntriesToFitnessLinkage.py T30_query_list linkage_file fitness_file
#                                                                            
# INPUTS:                                                                    
#  T30_query_list This is the list of entries to duplicate with their T30 tags
#  linkage_file && fitness_file, input files
#                                                                            
# OUTPUTS:
#  Script will rewrite new files, alongside the input files named:
#                _T30append_(date).txt
#                                                                            
#############################################################################
"""
    print HELP_TEXT
    return



################ MAIN FUNCTION
# imports and constants
import sys
import os
import datetime

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) != 4:
    print 'wrong number of arguments (try "-h" for help)'
    exit()

T30_query_list = sys.argv[1]
linkage_file = sys.argv[2]
fitness_file = sys.argv[3]
today = datetime.date.today()
append_str = '_T30append_' + str(today.year%100) + str(today.month) + str(today.day) + '.txt'

# Check input files
if not os.path.exists(T30_query_list):
    print 'T30_query_list "' + T30_query_list + '" does not exist'
    exit()
if not os.path.exists(linkage_file):
    print 'linkage_file "' + linkage_file+ '" does not exist'
    exit()
if not os.path.exists(fitness_file):
    print 'fitness_file "' + fitness_file + '" does not exist'
    exit()

# Construct and check output files
ext = linkage_file[-4:]
if not ext == '.txt':
    print 'linkage_file did\'t end in .txt, renaming will be ugly'
    linkage_output = linkage_file + append_str
else:
    linkage_output = linkage_file[:-4] + append_str

ext = fitness_file[-4:]
if not ext == '.txt':
    print 'fitness_file did\'t end in .txt, renaming will be ugly'
    fitness_output = fitness_file + append_str
else:
    fitness_output = fitness_file[:-4] + append_str

try:
    linkage_output_fid = open(linkage_output, 'w')
except:
    print 'Error opening linkage_output file: ' + linkage_output
    exit() 

try:
    fitness_output_fid = open(fitness_output, 'w')
except:
    print 'Error opening fitness_output file: ' + fitness_output
    exit() 

## Step 1: Hash the T30 queries so we know what we're looking for as we
#          scan through the other files.
T30 = set()
T30_count = 0
T30_fid = open(T30_query_list, 'r')
for line in T30_fid:
    T30_count += 1
    line = line.strip()
    assert(line[-3:] == 'T30') # check for T30
    line = line[:-3]              # and remove
    if line[-1] == '_':        # cases like YBL061C_T30, remove '_'
        line = line[:-1]
    T30.add(line)
T30_fid.close()

## Step 2: Scan through linkage file, output anything that doesn't match our T30 set
#  replicate anything that does with T30 on the end of the first column
linkage_count = 0
linkage_file_fid = open(linkage_file, 'r')
for line in linkage_file_fid:
    parsed_line = line.split('\t')
    if(parsed_line[0] in T30):
        linkage_count += 1
        if parsed_line[0][-1] == 'W' or parsed_line[0][-1] == 'C':
            parsed_line[0] = parsed_line[0]+'_T30'
        else:   # damp or tsq, no _ needed
            parsed_line[0] = parsed_line[0]+'T30'
        linkage_output_fid.write('\t'.join(parsed_line)) # write duplicate line

    linkage_output_fid.write(line) # preserve original line


linkage_output_fid.close()
linkage_file_fid.close()

## Step 3: Do the same for the fitness file
fitness_count = 0
fitness_file_fid = open(fitness_file, 'r')
for line in fitness_file_fid:
    parsed_line = line.split('\t')
    if(parsed_line[0] in T30):
        fitness_count += 1
        if parsed_line[0][-1] == 'W' or parsed_line[0][-1] == 'C':
            parsed_line[0] = parsed_line[0]+'_T30'
        else:
            parsed_line[0] = parsed_line[0]+'T30'
        fitness_output_fid.write('\t'.join(parsed_line)) # write duplicate line

    fitness_output_fid.write(line) # preserve original line


fitness_output_fid.close()
fitness_file_fid.close()

## Print a quick summary
print 'T30 input '+ str(T30_count) + '\n'
print 'T30 linkage '+ str(linkage_count) + '\n'
print 'T30 fitness '+ str(fitness_count) + '\n'

