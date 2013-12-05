import sys
import csv
import numpy
import getopt
import math
import random
import copy

def usage():
    print '$> python generatednadata.py <required args> [optional args]\n' + \
        '\t-c <#>\t\tNumber of clusters to generate\n' + \
        '\t-d <#>\t\tNumber of DNA strands per cluster\n' + \
        '\t-l <#>\t\tLength of DNA strands\n' + \
        '\t-o <file>\tFilename for the output of the DNA data\n'

def simDNA(d1, d2):
    '''
    Computes the similarity between two DNA strands.
    '''
    dlen=len(d1)
    sim = dlen
    for i in range(dlen):
        if d1[i] != d2[i] :
            sim = sim-1
    return sim

def tooSimilar(d, ds, maxSim):
    '''
    Computes the similarity between the DNA strand and all DNA strands
    in the list, and if any DNA strand in the list are closer than maxSim,
    this method returns true.
    '''
    for dna in ds:
        if simDNA(d, dna) > maxSim:
                return True
    return False

def handleArgs(args):
    # set up return values
    numClusters = -1
    numDNA = -1
    lenDNA = -1
    output = None

    try:
        optlist, args = getopt.getopt(args[1:], 'c:d:l:o:')
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    for key, val in optlist:
        if   key == '-c':
            numClusters = int(val)
        elif key == '-d':
            numDNA = int(val)
        elif key == '-l':
            lenDNA = int(val)
        elif key == '-o':
            output = val

    # check required arguments were inputted  
    if numClusters < 0 or numDNA < 0 or \
            lenDNA < 0 or \
            output is None:
        usage()
        sys.exit()
    return (numClusters, numDNA, lenDNA, output)

def drawOrigin(dlen):
    strand=[]
    bases=['A','C','G','T']
    for i in range(dlen):
        strand.append(bases[random.randint(0,3)])
    return strand



# start by reading the command line
numClusters, \
numDNA, \
lenDNA, \
output = handleArgs(sys.argv)

writer = csv.writer(open(output, 'w'))

# step 1: generate each centroid
centroids_radii = []
maxSim = 0.3*lenDNA
for i in range(0, numClusters):
    centroid_radius = drawOrigin(lenDNA)
    # is it far enough from the others?
    while (tooSimilar(centroid_radius, centroids_radii, maxSim)):
        centroid_radius = drawOrigin(lenDNA)
    centroids_radii.append(centroid_radius)

# step 2: generate the points for each centroid
minClusterDev = 0.1*lenDNA
maxClusterDev = 0.2*lenDNA
for i in range(0, numClusters):
    # compute the deviation for this cluster
    deviation = numpy.random.uniform(minClusterDev, maxClusterDev)
    cluster = centroids_radii[i]
    for j in range(0, numDNA):
        acentroid=copy.deepcopy(cluster)
        # generate number of differences
        numDiff = int(abs(numpy.random.normal(0, deviation)))
        # randomly choose locations to be altered
        indDiff=numpy.random.randint(0, lenDNA, numDiff)
        for k in range(0,numDiff):
            bases=['A','C','G','T']
            origBase=acentroid[indDiff[k]]
            bases.remove(origBase)
            # randomly alter location
            acentroid[indDiff[k]]=bases[random.randint(0,2)]
        # write the new strand out
        writer.writerow(acentroid)
