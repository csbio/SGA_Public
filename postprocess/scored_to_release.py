#!/usr/bin/env python
def help():
    HELP_TEXT = """
#############################################################################
# scored_to_relsease.py
#
# This script allows me to bypass the usual SGA post-processing pipeline
# and create a "realease" data file from a "scored" data file (E.G. CHI38)
#
# Author: Benjamin VanderSluis (bvander@cs.umn.edu)
# Revision: October 24, 2012
#
# USAGE:
# scored_to_release.py label scored_file > output
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
import fileinput

## Step 0: argument parsing and existence checks and help
if sys.argv.count('-h') + sys.argv.count('-help') + sys.argv.count('--help') > 0:
    help()
    exit()

if len(sys.argv) < 2:
    print 'too few arguments (try "-h" for help)'
    exit()

label       = sys.argv[1]
scored_file = sys.argv[2]

common_names = '/project/csbio/benjamin/Data/Master_Common_Ref_SGD.txt'
sys.stderr.write('using for common names:\n' + common_names + '\n')


common_fid = open(common_names, 'r')
cnames = {}
for line in common_fid:
	line = line.strip().split('\t')
	if len(line) > 1:
		cnames[line[0]] = line[1]
common_fid.close()


# Now ensure that these all exist and we're allowed to write the output
# if we fail because of this, we want to fail before doing a lot of work
if scored_file[-3:] == '.gz':
    scored_fid = fileinput.hook_compressed(scored_file, 'r')
else:
    scored_fid = open(scored_file, 'r')

header = ['Query ORF', 'Query gene name', 'Array ORF', 'Array gene name', 'Double mutant fitness', 'Epsilon', 'Standard deviation', 'P-value', 'Experiment']
print('\t'.join(header))

for line in scored_fid:
	line = line.strip().split('\t')
	newline = []

	# columns 1-4
	for i in range(2):
		orf = line[i]
		newline.append(orf)
		if '_' in orf:
			suff = orf[orf.find('_'):]
			orf = orf[0:orf.find('_')]
			if orf in cnames:
				newline.append(cnames[orf] + suff)
			else:
				newline.append(orf+suff)
	newline.append(line[10]) # DMF
	newline.append(str(float(line[10]) - float(line[9]))) # EPS
	newline.append(line[11]) # STD
	newline.append(line[4])
	newline.append(label)

	print('\t'.join(newline))

scored_fid.close()
