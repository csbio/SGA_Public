#!/usr/bin/env python3
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
# Revision: November 04, 2014                                                    
#                                                                            
# USAGE:                                                                     
# AddDoubleQueryFitnessToSMfile.py 
#                                                                            
# INPUTS:                                                                    
# ... <Parameters Hardcoded> ...
#  DM_query_file File with list of DM queries to add to the standard
#  DMF_file  Current DMF fitness standard. These values taken first
#  SM_input  This is the current SM standard, will be duplicated in output for a complete standard
#  SGA1 File with DM scores and deviations(output of compute_sgascore.m) (taken second)
#  SGA2, checked for dms if not found previously, list any number of such files...
#
# SGA_file can be gzipped (recommended)                                     
#
# OUTPUT: final query-fitness standard is written to stdout
#                                                                            
# NOTE: if ORF+YDL227C or YDL227C+ORF cannot be found in sga, we'll try to
#       replicate the ORF smf entry as the dmf entry
#
#############################################################################
"""
    print(HELP_TEXT)
    return

################ MAIN FUNCTION
# imports and constants
import sys
import os
import pandas as pd
import numpy as np
import ipdb

SGA_QUERY_COL = 0
SGA_ARRAY_COL = 1
SGA_DM_FIT_COL = 10
SGA_DM_STD_COL = 11
SGA_COLS = [SGA_QUERY_COL, SGA_ARRAY_COL, SGA_DM_FIT_COL, SGA_DM_STD_COL]
SGA_HEADER = ['query', 'array', 'dmf', 'dmf_std']

# Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    sys.exit()

if len(sys.argv) > 1:
    print('too many arguments (paths are hard-coded)')
    sys.exit()

DM_query_file = '/project/csbio/lab_share/SGA/rawdata/triple/150226/dm_queries.txt'
DMF_component_file = '/project/csbio/lab_share/SGA/Main/refdata/double_mutant_allele_composition_150227.csv'
DMF_std_file = '/project/csbio/benjamin/Triples/dmf/dmf_standard_150204.csv'
SM_input = '/project/csbio/lab_share/SGA/refdata/smf_t26_130417.txt'
SGA_files = ['/project/csbio/lab_share/SGA/Main/scored/131130/scored_sga_fg_t26_131130_scored_140103.txt',
'/project/csbio/lab_share/SGA/Main/scored/131130/scored_sga_ts_t26_131130_scored_140103.txt',
'/project/csbio/lab_share/SGA/Main/scored/131130/scored_sga_fg_t30_131130_scored_140103.txt',
'/project/csbio/lab_share/SGA/Main/scored/131130/scored_sga_ts_t30_131130_scored_140103.txt']

QUERY_FIT_OUTPUT = '/project/csbio/lab_share/SGA/refdata/smf_t26_130417_tr_150226.txt'

# Now ensure that we are not going to clobber the output
if os.path.exists(QUERY_FIT_OUTPUT):
    print('cowardly refusing to clobber outputfile\n')
    sys.exit()

## Step 1: build a dataframe to hold the DMF values #########################################################
# beginning with the list we need values for
DMF = pd.read_table(DM_query_file, header=None, names=['strain'])
DMF['orf'] = [x.split('_')[0] for x in DMF['strain']] 
DMF['tag'] = [x.split('_')[1] for x in DMF['strain']] 

## Step 2: load the DMF standard file and fill in matching values ############################################
DMF_experiment = pd.read_table(DMF_std_file, header=None, names=['strain', 'dmf', 'dmf_std'])
DMF_experiment = DMF_experiment.dropna()
DMF = pd.merge(DMF, DMF_experiment, how='left', on=['strain'])

## Step 3: Load the component table to prepare for SGA lookups
CT = pd.read_table(DMF_component_file)


## Step 4:
## Use the strain equiv (CT) to steal corresponding single mutants (at 26) as well
smf_dat = pd.read_table(SM_input, header=None, names=['strain', 'smf', 'smf_std'])
smf_dat['tag'] = [x.split('_')[1] for x in smf_dat['strain']]
smf_dat.index = smf_dat['tag']

missing = sum(np.isnan(DMF['dmf']))
print('initial missing ')
print(missing)
for i in np.where(np.isnan(DMF['dmf']))[0]:
    # if this is a YDL227C strain, we should find _tm in one of the SMstrain cols
    sm1_ix = CT['SM1strainID'].values == DMF['tag'][i]
    t = None
    if sm1_ix.any():
        t = CT['StrainID1'][sm1_ix].tolist()[0]
    else:
        # maybe no values (broadcasting?)
        sm2_ix = CT['SM2strainID'].values == DMF['tag'][i]
        if sm2_ix.any():
            t = CT['StrainID2'][sm2_ix].tolist()[0]

    if t != None:
        smf = smf_dat['smf'][smf_dat['tag'] == t]
        if len(smf) > 0:
            DMF['dmf'][i] = smf.values[0]
            DMF['dmf_std'][i] = smf.values[0]

missing = sum(np.isnan(DMF['dmf']))
print('post smf missing ')
print(missing)

# ipdb.set_trace() # is DMF set up right?

## Step 5: Iterate through the SGA_files and fill in any missing values ###################################
for sga_file in SGA_files:

    if sga_file[-2:] == '.gz':
        comp = 'gzip'
    else:
        comp = None
    sga_dat = pd.read_table(sga_file, compression=comp, usecols=SGA_COLS, names=SGA_HEADER)

    # create stripped orf versions
    sga_dat['q_tag'] = [x.split('_')[1] for x in sga_dat['query']]
    sga_dat['a_tag'] = [x.split('_')[1] for x in sga_dat['array']]

    # Go through missing values, convert _tm to _q _a pairs (not bothering with the reverse yet)
    # and check for this occurance in the current SGA file

    for i in np.where(np.isnan(DMF['dmf']))[0]:
        # check if its a real DM strain
        if DMF['tag'][i] in CT['DMstrainID'].values:
            # get the corresponding tags, and pull dmf from SGA
            ct_ix = CT['DMstrainID'] == DMF['tag'][i]
            q = CT['StrainID1'][ct_ix].values[0] 
            a = CT['StrainID2'][ct_ix].values[0]

            # search one at a time as a sort of short-circuit
            sga_q = sga_dat[sga_dat['q_tag'] == q]
            sga_qa = sga_q[sga_q['a_tag'] == a]

            if len(sga_qa) > 1:
                pdb.set_trace()

            if len(sga_qa) > 0:
                DMF['dmf'][i] = np.nanmean(sga_qa['dmf'])
                DMF['dmf_std'][i] = np.nanmean(sga_qa['dmf_std'])

## Step 6: Print out the result
missing = sum(np.isnan(DMF['dmf']))
print('post sga missing ')
print(missing)
# at home I had to use "columns" and on site "cols"
# DMF.to_csv(QUERY_FIT_OUTPUT, na_rep='NaN', sep='\t', index=False, columns=['strain', 'dmf', 'dmf_std'])
DMF.to_csv(QUERY_FIT_OUTPUT, na_rep='NaN', sep='\t', header=False, index=False, cols=['strain', 'dmf', 'dmf_std'])
print('finished')

## Step 7: cat the smf file to the DMF file 
execl = 'cat ' + SM_input + ' >> ' + QUERY_FIT_OUTPUT
os.system(execl)
