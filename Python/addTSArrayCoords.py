#ek
# Makes a new coordinate file by splitting gene names on "_" and matching
#   the name before the underscore to the gene names in the original 
#   coordinates file.
# Created: 9/7/11

# Filenames
newNamesSgaFilename = '/project/csbio/lab_share/SGA/rawdata'+\
                      '/TS-Array-v6_dump_110722-13-12-36.txt'
newCoordsFilename = '/project/csbio/lab_share/SGA/refdata/'+\
                    'chrom_coordinates_wTS_110907.txt'

# Get all arrays from raw data
allArrays = {}
with newNamesSgaFile as f:
    for line in f:
        allArrays[line.split()[1]] = 1
allArrays = allArrays.keys()

# Get coordinate data for all orfs
coordsDict = {}
with open('/project/csbio/lab_share/SGA/refdata/chrom_coordinates.txt') as f:
    coordLines = f.readlines()
    for line in coordLines:
        L = line.split('\t', 1)
        coordsDict[L[0]] = L[1]

# Make dictionary with coordinate data for each ts array gene
tsCoordsDict = {}
for arrayAllele in allArrays:
    arraySplit = arrayAllele.split('_', 1)
    if len(arraySplit) == 2:
        try:
            tsCoordsDict[arrayAllele] = coordsDict[arraySplit[0]]
        except KeyError:
            print 'No coord data for:', arrayAllele

# Write file with all coords data
newCoordLines = coordLines + ['%s\t%s' % x for x in tsCoordsDict.items()]
with open newCoordsFilename as fout:
    fout.writelines(newCoordLines)
