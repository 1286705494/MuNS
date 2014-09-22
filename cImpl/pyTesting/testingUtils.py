
import sys
import csv
import math
import numpy as np
from operator import itemgetter



def loadPart(sAbsPartFile, sMode):
    """
    Load the partitions from file.
    """
    
    llPartitions = []
    
    if sMode == 'partition':
        with open(sAbsPartFile, 'r') as fPart:
            # read in the actual partition file
            csvPart = csv.reader(fPart, delimiter=',')
            for lRow in csvPart:
                llPartitions.append([int(x) for x in lRow])
    elif sMode == 'membership':
        currVert = 0
        with open(sAbsPartFile, 'r') as fPart:
            csvPart = csv.reader(fPart, delimiter=',')
            for lRow in csvPart:
                # only take first element as partition index
                pos = int(lRow[0])
                # add partitions while we still don't have enough
                while len(llPartitions) <= pos:
                    llPartitions.append([]) 
                llPartitions[pos].append(currVert)
                
                currVert += 1            
    elif sMode == 'list':
        currVert = 0
        
        with open(sAbsPartFile, 'r') as fPart:
            csvPart = csv.reader(fPart, delimiter=',')
            for lRow in csvPart:
                # only take first element as partition index
                vertId = int(lRow[0]) # we dont' use
                pos = int(lRow[1])
                # add partitions while we still don't have enough (assume pos start at 0)
                while len(llPartitions) <= pos:
                    llPartitions.append([]) 
                llPartitions[pos].append(currVert)
                
                currVert += 1        
            
    return llPartitions


############################################################################################

def avgIntraSim(llSim, llPartitions):
    """
    Compute the average intra-partition similarity for the specified partitions.
    """
    
    totalSim = 0
    
    totalPairs = 0
    
    for lPart in llPartitions:
        # calculate the number of element
        totalPairs += len(lPart) * len(lPart)
        for x in lPart:
            for y in lPart:
                totalSim += llSim[x][y]
                
    return totalSim / totalPairs


def avgInterSim(llSim, llPartitions):
    """
    Compute the average inter-partition similarity for the specified partitions.
    """    
    
    totalSim = 0
    
    totalPairs = 0
    
    partNum = len(llPartitions)
    for p in range(0, partNum):
        lPart1 = llPartitions[p]
        for j in range(p+1, partNum):
            lPart2 = llPartitions[j]
            totalPairs += len(lPart1) * len(lPart2)
            for x in lPart1:
                for y in lPart2:
                    totalSim += llSim[x][y]
                
    return (totalSim / totalPairs)


def avgRank(lSim, rowLen, llPartition, epsilon):
    """
    Compute the average similarity ranking within the specified partitions.
    
    @param epsilon: If two values in llSim < this, then it is considered as being the same.
    """
    
    
    # covert to a vector
#     lSim = [item for lSublist in llSim for item in lSublist]
#     rowLen = len(llSim[0])
#     print lSim
    
#     # sort the similarities with indices into it
#     lSortedIndices = [i for (i,val) in sorted(enumerate(lSim), key=itemgetter(1), reverse=False)]
# #     print lSortedIndices
#     
#     lRank = [0 for x in range(0,len(lSortedIndices))]
#     
#     # compute the rank
#     for (i,ind) in enumerate(lSortedIndices):
#         lRank[ind] = float(i+1)  / len(lRank)
# #     print lRank

    lRank = percentileRank(lSim, epsilon)
        
    
    totalIntraPairs = 0
    totalIntraRank = 0;
    for lPart in llPartition:
        # calculate the number of element
        totalIntraPairs += len(lPart) * len(lPart)
        for x in lPart:
            for y in lPart:
                totalIntraRank += lRank[x*rowLen + y]
                
                
    totalInterRank = 0
    totalInterPairs = 0
    
    partNum = len(llPartition)
    for p in range(0, partNum):
        lPart1 = llPartition[p]
        for j in range(p+1, partNum):
            lPart2 = llPartition[j]
            totalInterPairs += 2*len(lPart1) * len(lPart2)
            for x in lPart1:
                for y in lPart2:
                    totalInterRank += lRank[x*rowLen + y]
                    totalInterRank += lRank[y*rowLen + x]
                
#     print totalIntraRank
#     print totalInterRank
#     print totalIntraPairs
#     print totalInterPairs
    return float(totalIntraRank) / totalIntraPairs, float(totalInterRank) / totalInterPairs               


def percentileRank(lSim, epsilon):
    """
    Compute the percentile rank of llSim.
    
    @param epsilon: If two values in llSim < this, then it is considered as being the same.
    """
    
    # compute histogram
    binNum = 1.0 / epsilon
#     binBoundaries = range(0, binNum)
    
    lHist, binBoundaries = np.histogram(lSim, bins=binNum)
#     print list(lHist)
#     print list(binBoundaries)
    vHistCumsum = np.cumsum(lHist)
    
#     print vHistCumsum
    
    itemNum = len(lSim)
    lRank = [0 for x in lSim]
    for (i, sim) in enumerate(lSim):
        bin = findBin(sim, binBoundaries)
#         print str(sim) + " -> " + str(bin)
        if bin <= 1:
            lRank[i] = lHist[bin-1] / (2.0 * itemNum)
        else:
            lRank[i] = (vHistCumsum[bin-2] + lHist[bin-1] / 2.0) / itemNum 
        
    return lRank
    
    
    
    
     
                
    


def avgPrecisionAtK(llSim, llPartition, topN):
    """
    Compute the precision of the topN peers of each vertex and compare with the vertices in the same class as each vertex.
    """
    
    totalPrecision = 0
    
    elemNum = 0
    
    # loop through the elements in llPartition
    for lPart in llPartition:
        setPart = set(lPart)
        for v in lPart:
            # get corresonding row for sim
            lSortedIndices = [i for (i,val) in sorted(enumerate(llSim[v]), key=itemgetter(1), reverse=True)]
            setTopN = set(lSortedIndices[:topN])
                   
            intersection = len(setTopN.intersection(setPart))
            totalPrecision += float(intersection) / len(setTopN)
            
            elemNum += 1
             
    return totalPrecision / elemNum 
    
    
        
def avgRPrecision(llSim, llPartition):
    """
    Compute the precision of the topN peers of each vertex and compare with the vertices in the same class as each vertex.
    """
    
    totalPrecision = 0
    
    elemNum = 0
    
    
    # loop through the elements in llPartition
    for lPart in llPartition:
        setPart = set(lPart)
        topN = len(setPart)
        for v in lPart:
            # get corresonding row for sim
            lSortedIndices = [i for (i,val) in sorted(enumerate(llSim[v]), key=itemgetter(1), reverse=True)]
            setTopN = set(lSortedIndices[:topN])
                   
            intersection = len(setTopN.intersection(setPart))
            totalPrecision += float(intersection) / topN
            
            elemNum += 1
             
    return totalPrecision / elemNum

    
def avgRetrievalStats(llSim, llPartition, topN):
    """
    Compute the precision of the topN peers of each vertex and compare with the vertices in the same class as each vertex.
    """
    
    totalPrecision = 0
    totalRecall = 0
    totalFscore = 0
    
    elemNum = 0
    
    # loop through the elements in llPartition
    for lPart in llPartition:
        setPart = set(lPart)
        for v in lPart:
            # get corresonding row for sim
            lSortedIndices = [i for (i,val) in sorted(enumerate(llSim[v]), key=itemgetter(1), reverse=True)]
            setTopN = set(lSortedIndices[:topN])
                   
            intersection = len(setTopN.intersection(setPart))
            precision = float(intersection) / len(setTopN)
            recall = float(intersection) / topN
            fscore = 0
            if precision > 0 and recall > 0:
                fscore = float(2 * precision * recall) / (precision + recall)
            
            totalPrecision += precision
            totalRecall += recall
            totalFscore += fscore             
        
            elemNum += 1
             
    return totalPrecision / elemNum, totalRecall / elemNum, totalFscore / elemNum


def avgMAP(llSim, llPartition):
    """ 
    Mean average precision.
    """
    
    totalPrecision = 0

    
    elemNum = 0
    
    # loop through the elements in llPartition
    for lPart in llPartition:
        setPart = set(lPart)
        
        for v in lPart:
            currPrecision = 0
            # get corresonding row for sim
            lSortedIndices = [i for (i,val) in sorted(enumerate(llSim[v]), key=itemgetter(1), reverse=True)]
            
            # for each vertex in the same partition as v (apart from v), compute the precision of retrieving that element
            for o in lPart:
                # find its position 
                if o != v:
                    pos = lSortedIndices.index(o)
                    setTopN = set(lSortedIndices[:pos])
                   
                    intersection = len(setTopN.intersection(setPart))
                    precision = 0
                    if len(setTopN) > 0:
                        precision = float(intersection) / len(setTopN)
                    assert(precision <= 1)
            
                    currPrecision += precision
            
            elemNum += 1
            
            if len(lPart)-1 > 0:
                avgPartPrecision = float(currPrecision) / (len(lPart)-1)
                assert(avgPartPrecision <= 1)
                totalPrecision += avgPartPrecision
                
        
             
    return totalPrecision / elemNum      


#############################################################


def computeCompactness(mDis, vClusters, vMedriods):
    
    compactness = 0
    
    totalDis = 0
    for (v, clus) in enumerate(vClusters):
        totalDis += mDis[v,clus]
        
    totalMedriodDis = 0
    for i in range(0, len(vMedriods)):
        for j in range(i+1, len(vMedriods)):
            totalMedriodDis += mDis[vMedriods[i], vMedriods[j]]
    
    compactness = float(totalDis) / totalMedriodDis
    
    return compactness


#################################################################

def computeDShell(mDiGraph, outDegBinNum, inDegBinNum):
    """
    Compute the Directed-Shell of mDiGraph.
    maxK - maximum outdegree value 
    maxL - maximum indegree value
    stepOutDeg - step size per bin for out degree.
    stepInDeg - step size per bin for in degree.
    """
    
    # compute the in and out degree of each vertex
    hOutDeg = mDiGraph.out_degree()
    hInDeg = mDiGraph.in_degree()
    
#     print hOutDeg
#     print hInDeg
    
    # construct bins                     
    lHist, lOutBinEdges = np.histogram(np.array(hOutDeg.values()), bins=outDegBinNum)
    lHist, lInBinEdges = np.histogram(np.array(hInDeg.values()), bins=inDegBinNum)
#     print lOutBinEdges
#     print lInBinEdges
    
    currInBinIndex = 0
    currOutBinIndex = 0
    
    llDShells = [[] for x in range(outDegBinNum * inDegBinNum)]
    
    # construct the D-Shells
    for (vertIndex, outDeg) in hOutDeg.items():
        inDeg = hInDeg[vertIndex]
        
        outBin = findBin(outDeg, lOutBinEdges)
        inBin = findBin(inDeg, lInBinEdges)
#         print outBin
#         print inBin
            
        llDShells[(outBin-1) * inDegBinNum + (inBin-1)].append(vertIndex)
        
    return llDShells
        
    
        
        
def findBin(deg, lBinEdges):
    
    bin = 1
    while bin <= len(lBinEdges):
        if deg <= lBinEdges[bin]:
            break
        bin += 1
        
    return bin
    
    
def intraCrossShellsRank(llDShells, lRank, rowNum, colNum):
    """
    Compute the intra and cross shell average rank.
    """
    
    mIntraRank = np.zeros((rowNum, colNum))
    
    # intra for each dshell
    for r in range(rowNum):
        for c in range(colNum):
            lPart = llDShells[r * colNum + c]
    
            totalIntraRank = 0;
            partSize = len(lPart)
            if partSize > 0:
                for i in range(0, partSize):
                    for j in range(i+1, partSize):
                        totalIntraRank += lRank[lPart[i]*colNum + lPart[j]]
                        totalIntraRank += lRank[lPart[j]*colNum + lPart[i]]    
    
                # diagonal
                for x in lPart:
                    totalIntraRank += lRank[x * colNum + x]
                
                # change if we remove diagonal
                mIntraRank[r,c] = totalIntraRank / float(partSize * partSize)
            
    
    
    # shell difference
    mShellInterRank = np.zeros((rowNum, colNum))
    mShellInterPairNum = np.zeros((rowNum, colNum))

    for (i, lPartOuter) in enumerate(llDShells):
        if len(lPartOuter) > 0:
            for j in range(i+1, len(llDShells)):
                lPartInner = llDShells[j]
                if len(lPartInner) > 0:
                    # get indicues
                    outerR = math.floor(i / colNum)
                    outerC = i % colNum
                    
                    innerR = math.floor(j / colNum)
                    innerC = j % colNum
                    
                    mShellInterPairNum[abs(innerR - outerR), abs(innerC - outerC)] += len(lPartOuter) * len(lPartInner)          
                    for x in lPartOuter:
                        for y in lPartInner:
                            mShellInterRank[abs(innerR - outerR), abs(innerC - outerC)] += lRank[x * colNum + y] 

    
    # do row first
#     for outerR in range(rowNum):
#         for outerC in range(colNum):
#             print 'outer = ' + str(outerR) + ',' + str(outerC)
#             lPartOuter = llDShells[outerR * colNum + outerC]
#             if (len(lPartOuter) == 0):
#                 continue
#             print lPartOuter
#             
#             # loop through the other shells
#             for innerR in range(outerR, rowNum):
#                 for innerC in range(outerC, colNum):
#                     
#                     if innerR == outerR and innerC == outerC:
#                         continue
#                     print 'inner = ' + str(innerR) + ',' + str(innerC)
#                     
#                     
#                     lPartInner = llDShells[innerR * colNum + innerC]
#                     if (len(lPartInner) == 0):
#                         continue;
#                     print lPartInner
#                     print innerR - outerR
#                     print innerC - outerC
#                     mShellInterPairNum[innerR - outerR, innerC - outerC] += len(lPartOuter) * len(lPartInner) 
#                     for x in lPartOuter:
#                         for y in lPartInner:
#                             mShellInterRank[innerR - outerR, innerC - outerC] += lRank[x * colNum + y] 
                            
                            
            
    
    # normalise each entry of mShellInterRank
    for r in range(rowNum):
        for c in range(colNum):
            if mShellInterPairNum[r,c] > 0:
                mShellInterRank[r,c] /= float(mShellInterPairNum[r,c])
    
    return mIntraRank, mShellInterRank
    
    
    