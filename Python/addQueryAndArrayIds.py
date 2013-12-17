#!/usr/bin/env python
def help():
    HELP_TEXT = """ 
#############################################################################
# addQueryAndArrayIds.py
# Using for appending annotations to gene names for smf collections
#
# Author: Elizabeth Koch (koch@cs.umn.edu)
# Revision: November 29, 2011
#
#                    Now, assumes inFile is tab-delimited and
#                    changes only the first two fields:
#
#   "query    array   rest" --> "query_prefixId array_prefixId rest"
#
# USAGE: addQueryAndArrayIds.py ...
# query_id_file array_id_file input_file output_file query_prefix array_prefix
# refdata/sn_Ids_111129.csv refdata/dma_Ids_111129.csv input_file output_file sn dma
#
#############################################################################
"""
    print HELP_TEXT
    return

import sys

def getIdDict(idFilename):
    orfToId = {}
    with open(idFilename) as f:
        for line in f:
            L = line.split('\t')
            if L[4] not in orfToId:
                orfToId[L[4]] = L[1]
    return orfToId

if __name__ == '__main__':
    if len(sys.argv) != 7:
        help()
        sys.exit(0)
    else:
        idFilename_q = sys.argv[1]
        idFilename_a = sys.argv[2]
        inFilename = sys.argv[3]
        outFilename = sys.argv[4]
        prefix_q = sys.argv[5]
        prefix_a = sys.argv[6]

    orfToId_q = getIdDict(idFilename_q)
    orfToId_a = getIdDict(idFilename_a)

    fin = open(inFilename)
    fout = open(outFilename, 'w')
    for line in fin:
        (query, array, rest) = line.split('\t', 2)

        if '_' not in query:
            if query in orfToId_q:
                query = query + '_' + prefix_q + orfToId_q[query]

        if '_' not in array:
            if array in orfToId_a:
                array = array + '_' + prefix_a + orfToId_a[array]

        fout.write('%s\t%s\t%s' % (query, array, rest))

    fin.close()
    fout.close()
