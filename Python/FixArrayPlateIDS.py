#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# FixArrayPlateIDS.py
#
# This script tranforms array plate ids
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: June 30, 2011
#
# USAGE:
# FixArrayPlateIDS.py input_file output_file
#
# 357 -> 1
# 358 -> 2
# 359 -> 3
# 360 -> 4
# 361 -> 5
# 362 -> 6
# 363 -> 7
# 364 -> 8
# 365 -> 9
# 366 -> 10
# 367 -> 11
# 368 -> 12
# 369 -> 13
# 370 -> 14
#
# INPUTS:
#
# OUTPUTS:
#
#############################################################################
"""
    print HELP_TEXT
    return

################ MAIN FUNCTION
# imports and constants
import sys
import os
# import fileinput

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) < 3:
    print 'too few arguments (try "-h" for help)'
    exit()

input_file = sys.argv[1]
output_file = sys.argv[2]

# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if not os.path.exists(input_file):
    print 'input_file "' + input_file + '" does not exist'
    exit()

try:
    output_fid = open(output_file, 'w')
except:
    print 'Error opening output file ' + output_file
    exit() 


# Step one, the hash
Valid_set = set(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14'])
ID_HASH = { '357':'1', '358':'2', '359':'3', '360':'4', '361':'5', '362':'6', '363':'7', '364':'8', '365':'9', '366':'10', '367':'11', '368':'12', '369':'13', '370':'14'}

input_fid = open(input_file, 'r')
for line in input_fid:
    line = line.split('\t')
    if line[2] in Valid_set:
        output_fid.write('\t'.join(line))
    elif line[2] in ID_HASH:
        line[2] = ID_HASH[line[2]]
        output_fid.write('\t'.join(line))
    else:
        exit('unknown plate id found')

input_fid.close()
output_fid.close()
