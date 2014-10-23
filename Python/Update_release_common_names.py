#!/usr/bin/env python3
def help():
    HELP_TEXT = """
#############################################################################
#  Update_release_common_names.py
#
# This script updates columns 2,4 in a "release" file using
# the provided common name mapping. All Common names coming
# out will come from the file. So if there's one in the input
# which no longer exists, it will get lost (on purpose).
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: June 25, 2014
#
# USAGE:
# Update_release_common_names.py inputfile common_map > outputfile
#
# INPUTS:
#    inputfile: and sga file in "release" format
#    common_map: <ORF>\t<COMMON> format
#     e.g. /Data/YeastGeneMap/Orf_Common_Map_140617.txt
#
# OUTPUTS:
#     STD_OUT
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

if len(sys.argv) < 3:
    print('too few arguments (try "-h" for help)')
    exit()

inputfile = sys.argv[1]
cmapfile = sys.argv[2]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(inputfile):
    print('inputfile "' + inputfile + '" does not exist')
    exit()
if not os.path.exists(cmapfile):
    print('cmapfile "' + cmapfile + '" does not exist')
    exit()

cmapfid = open(cmapfile, 'r')
inputfid = fileinput.hook_compressed(inputfile, 'r')

# Step 1, load the cmap into a dict
cmap = {}
for line in cmapfid:
    line = line.strip().split('\t')
    if len(line) > 1 and len(line[1]) > 0:
        cmap[line[0]] = line[1]

cmapfid.close()
# Step 2, rip through the file, replace cols 2,4 and spit to STDout
# header should go through fine

for line in inputfid:
    line = line.strip().split('\t')
    query = line[0].split('_')[0]
    array = line[2].split('_')[0]
    if query in cmap:
        line[1] = cmap[query]
    else:
        line[1] = query

    if array in cmap:
        line[3] = cmap[array]
    else:
        line[3] = array

    line = '\t'.join(line)
    print(line),

inputfid.close()





