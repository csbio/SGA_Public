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
# AddDoubleQueryFitnessToSMfile.py DM_query_file DMF_component_file DMF_file SMF_file SGA1 SGA2 ...
#                                                                            
# INPUTS:                                                                    
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

## Step 0: argument parsing and existence checks and help
# if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
#     help()
#     sys.exit()

# if len(sys.argv) < 6:
#     print('too few arguments (try "-h" for help)')
#     sys.exit()

# DM_query_file = sys.argv[1]
# DMF_component_file = sys.argv[2]
# DMF_std_file = sys.argv[3]
# SM_input = sys.argv[4]
# SGA_files = sys.argv[5:]

# DM_query_file = '/project/csbio/lab_share/SGA/rawdata/triple/141021/dm_queries'
# DMF_component_file = '/project/csbio/lab_share/SGA/Main/refdata/double_mutant_allele_composition_141104.csv'
# DMF_std_file = '/project/csbio/lab_share/SGA/rawdata/triple/dmf/dmf_standard_141103.csv'
# SM_input = '/project/csbio/lab_share/SGA/refdata/smf_t26_130417.txt'
# SGA_files = ['/project/csbio/lab_share/SGA/Main/scored/131130/scored_sga_fg_t26_131130_scored_140103.txt',
# '/project/csbio/lab_share/SGA/Main/scored/131130/scored_sga_ts_t26_131130_scored_140103.txt',
# '/project/csbio/lab_share/SGA/Main/scored/131130/scored_sga_fg_t30_131130_scored_140103.txt',
# '/project/csbio/lab_share/SGA/Main/scored/131130/scored_sga_ts_t30_131130_scored_140103.txt']

DM_query_file = '/heap/data/stage/dmf/data/dm_queries'
DMF_component_file = '/heap/data/stage/dmf/data/double_mutant_allele_composition_141104.csv'
DMF_std_file = '/heap/data/stage/dmf/data/dmf_standard_141105.csv'
SM_input = '/heap/data/stage/dmf/data/smf_t26_130417.txt'
SGA_files = ['/heap/data/stage/dmf/data/scored_sga_fg_t26_131130_scored_140103.txt',
'/heap/data/stage/dmf/data/scored_sga_ts_t26_131130_scored_140103.txt',
'/heap/data/stage/dmf/data/scored_sga_fg_t30_131130_scored_140103.txt',
'/heap/data/stage/dmf/data/scored_sga_ts_t30_131130_scored_140103.txt']

# QUERY_FIT_OUTPUT = '/project/csbio/lab_share/SGA/refdata/smf_t26_130417_tm_141104.txt'
QUERY_FIT_OUTPUT = '/home/slice/stage/dmf/smf_t26_130417_tr_141105.txt'

# Now ensure that these all exist and we're allowed to write the output
for check_file in sys.argv[1:]:
    if not os.path.exists(check_file):
        print(check_file + ' does not exist', file=sys.stderr)
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


## I should use the strain equiv to steal single mutants (at 26) as well
## Step 4:
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
    # sga_dat['query_orf'] = [x.split('_')[0] for x in sga_dat['query']]
    # sga_dat['array_orf'] = [x.split('_')[0] for x in sga_dat['array']]
    sga_dat['q_tag'] = [x.split('_')[1] for x in sga_dat['query']]
    sga_dat['a_tag'] = [x.split('_')[1] for x in sga_dat['array']]

    # Go through missing values, convert _tm to _q _a pairs (not bothering with the reverse yet)
    # and check for this occurance in the current SGA file
    # there's a "pandas way" to do this, but time is a factor

    # for i,row_vals in DMF[np.isnan(DMF['dmf'])].iterrows():
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
DMF.to_csv(QUERY_FIT_OUTPUT, sep='\t', index=False, columns=['strain', 'dmf', 'dmf_std'])
print('finished')

## Step 7: cat the smf file to the DMF file 
execl = 'cat ' + SM_input + ' >> ' + QUERY_FIT_OUTPUT
os.system(execl)
