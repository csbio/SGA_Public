#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# AddDoubleQueryFitnessToSMfile.py                                           
#                                                                            
# This script gleans double mutant information from a scored SGA file and    
# adds this information to the SM standard file for the processing of       
# double-mutant queries.                                                    
#                                                                            
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)                          
# Revision: December 03, 2012                                                    
#                                                                            
# USAGE:                                                                     
# AddDoubleQueryFitnessToSMfile.py SM_input SM_output SGA_file [DM_query_file]
#                                                                            
# INPUTS:                                                                    
#  SM_input  This is the current SM standard, will be duplicated in output   
#  SM_output New standard with old data plus DM values for double queries    
#  SGA_file  File with DM scores and deviations (in columns xx and yy) (output of compute_sgascore.m)
#  DM_query_file File with list of DM queries in first column, if this       
#                 argument is ommitted, the SGA_file is scanned for queries  
#                 with a '+' in the name. (slower)                           
#                                                                            
# SGA_file and/or DM_query_file can be gzipped (recommended)                                     
#
# NOTE: if ORF+YDL227C or YDL227C+ORF cannot be found in sga, we'll try to
#       replicate the ORF smf entry as the dmf entry
#
# TODO: iff ORF1 and ORF2 are in the fitness file, but they haven't been
# screened against each other, output their product?
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

SGA_QUERY_COL = 0
SGA_ARRAY_COL = 1
SGA_DM_FIT_COL = 10
SGA_DM_STD_COL = 11

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) < 5:
    print 'too few arguments (try "-h" for help)'
    exit()

SM_input = sys.argv[1]
SM_output = sys.argv[2]
SGA_file = sys.argv[3]
if len(sys.argv) > 4:
    DM_query_file = sys.argv[4]
    QUERIES_IN_SGA = False
else:
    # use the SGA file to locate queries
    DM_query_file = sys.argv[3]
    QUERIES_IN_SGA = True

# Now ensure that these all exist and we're allowed to write the output
if not os.path.exists(SM_input):
    print 'SM_input file "' + SM_input + '" does not exist'
    exit()
if not os.path.exists(SGA_file):
    print 'SGA_file "' + SGA_file + '" does not exist'
    exit()
if not os.path.exists(DM_query_file):
    print 'DM_query_file "' + DM_query_file + '" does not exist'
    exit()
try:
    SM_output_fid = open(SM_output, 'w')
except:
    print 'Error opening SM_output file: ' + SM_output
    exit() 


## Step 1: assemble a list of DM queries for which we need fitness data
# if no seperate query file is given we'll do this later as we scan SGA
DM_Fitness = {}; 
if not QUERIES_IN_SGA:
    DM_query_fid = open(DM_query_file, 'r')
    for line in DM_query_fid:
        line = line.strip()
        # Remove any strain annotation on a DM query
        line = line.split('_')[0]
        DM_Fitness[line] = -1
    DM_query_fid.close()
else:
    DM_query_fid = fileinput.hook_compressed(DM_query_file, 'r')

    for line in DM_query_fid:
        line = line.split('\t')
        if line[SGA_QUERY_COL].find('+') > 0:
            DM_Fitness[line[SGA_QUERY_COL]] = -1
    DM_query_fid.close()


## Step 2: Scan through the SGA_file and update fitness values in the dict
SGA_fid = fileinput.hook_compressed(SGA_file, 'r')

# for now, entries in the fitness file are annotated
# but double mutant queries are not, so we must .split()
# to find them
for line in SGA_fid:
    line = line.strip().split()
    query = line[SGA_QUERY_COL].split('_')[0]
    array = line[SGA_ARRAY_COL].split('_')[0]

    AB = '+'.join([query, array])
    BA = '+'.join([array, query])

    if BA in DM_Fitness:
       AB = BA # switch the DM query name, appy same logic

    if AB in DM_Fitness: 
        if DM_Fitness[AB] == -1:
            # this is a DM query we haven't seen
            DM_Fitness[AB] = (line[SGA_DM_FIT_COL], line[SGA_DM_STD_COL])
        else:
            # We have seen AB (meaning we've already hit this pair as BA)
            fit = (float(DM_Fitness[AB][0]) + float(line[SGA_DM_FIT_COL])) / 2
            std = (float(DM_Fitness[AB][1]) + float(line[SGA_DM_STD_COL])) / 2
            DM_Fitness[AB] = (str(fit), str(std))
SGA_fid.close()


## Step 3: Spit out the original SM_fitness file, and append our new queries
# We may need some of these again so we'll hash them
# watch out for '_' ie ORF_sn...
SM_input_fid = open(SM_input, 'r')
smf_hsh = {}
for line in SM_input_fid:
    gene = line.strip().split()
    smf_hsh[gene[0].split('_')[0]] = '\t'.join(gene[1:])
    if gene[0] in DM_Fitness:
        # don't print lines we're going to print later
        continue 

    SM_output_fid.write(line)
SM_input_fid.close()

for key in DM_Fitness:
    if DM_Fitness[key] == -1:
        print key + ' notfound in SGA_file'
        single = ''
        if '+' not in key:
            error('+ not found in ' + key)

        [a, b] = key.split('+')
        if '_' in a:
            a = a[:a.find('_')]
        if '_' in b:
            b = b[:b.find('_')]

        if a == 'YDL227C':
            single = b
        elif b == 'YDL227C':
            single = a

        if single != '' and single in smf_hsh:
            print '   using '+single+' SMF for ' + key  
            SM_output_fid.write(key + '\t' + smf_hsh[single] + '\n')
        else:
            print '   '+single + ' notfound in SMF_file'

    else:    
        print key + ' found in SGA_file'
        (fit,std) = DM_Fitness[key]
        SM_output_fid.write('\t'.join([key, fit, std])+'\n')

SM_output_fid.close()
