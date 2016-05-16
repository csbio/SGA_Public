#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# prep_chi_queries.py
#
# This script removes trailing tags from the queries in a chi file.
# e.g. YGR245C_tsq523_michael_chi -> YGR245C_tsq523
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: June 30, 2014
#
# USAGE:
# prep_chi_queries.py inputfile > outputfile
#
# INPUTS: Gzipped chi scorefile
#
# OUTPUTS: STDOUT
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

inputfile = sys.argv[1]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(inputfile):
    print('inputfile "' + inputfile + '" does not exist')
    exit()

# compressed hook example, autodetects...
inputfid = fileinput.hook_compressed(inputfile, 'r')

for line in inputfid:
    line = line.strip().split('\t')
    query = line[0].split('_')
    line[0] = '_'.join(query[:2])
    print('\t'.join(line))
