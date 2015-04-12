'''
Created on 22/03/2015

@author: Jeffrey Chan, 2015
'''

import sys
import getopt
import glob
import csv
import os.path
import string
import subprocess as sp
import numpy as np
import networkx as nx
from testingUtils import *
from itertools import chain

def main():
    # process options
    try:
        # option list
        options = "p:t:e:b:n:l:s:i:a:d:hm:"
        # get options
        optList, remainArgs = getopt.gnu_getopt(sys.argv[1:], options)
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
        usage(sys.argv[0])
    
     
    sMode = 'list'
    maxIter = 20
    convEpsilon = 0.03
    ioBalance = 0.5
    topNCompared = 10
    simEpsilon = 0.001
    earlySimStopThres = 0.1
    icebergThres = 0.8
    iceBergAproxFactor = 0.5
    bShellAnalysis = False
    inBinNum = 10
    outBinNum = 10    
    for opt, arg in optList:
        if opt == "-p":
            sMode = arg
        elif opt == '-t':
            maxIter = int(arg)
        elif  opt == '-e':
            convEpsilon = float(arg)
        elif  opt == '-b':
            ioBalance = float(arg)
        elif  opt == '-n':
            topNCompared = int(arg)
        elif  opt == '-l':
            simEpsilon = float(arg)
            assert(simEpsilon <= 1.0 and simEpsilon >= 0.0)            
        elif opt == '-s':
            earlySimStopThres = float(arg)
        elif opt == '-i':
            icebergThres = float(arg)       
        elif opt == '-a':
            iceBergAproxFactor = float(arg)
        elif opt == "-h":
            bShellAnalysis = True            
        elif  opt == '-m':
            lBin = string.split(arg, ',')
            assert(len(lBin) == 2)
            outBinNum = int(lBin[0])
            inBinNum = int(lBin[1])               
        else:
            print >> sys.stderr, sys.argv[0] + " -m <mat index wanted> [edge list]"
            sys.exit(2)    
    
    # check number of arguments
    if len(remainArgs) < 4:
        usage(sys.argv[0])
    
    fileNum = int(remainArgs[0])
    lsGraphFile = remainArgs[1:fileNum+1]    
    lsPartFile = remainArgs[fileNum+1:2*fileNum+1]
    sResultsOutPrefix = remainArgs[2*fileNum+1]
#     vertNum = int(remainArgs[3])
    
    print lsGraphFile
    
#     lLambdaRange = [0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]
#     lLambdaRange = [0.8, 0.9, 0.95, 1.0]
    lLambdaRange = [0.90]
#     lBetaRange = [0.5]
#     lLambdaRange = [0.9]
#     sInitAlgor = 'binaryInit'
    sInitAlgor = 'degRatioInit'
#     sInitAlgor = 'degBinaryInit'
    ioBalance = 0.5
#     lsAlgor = ['autosim', 'rolesim', 'matchsim', 'simrank', 'prank']
    lsAlgor = ['autosim', 'rolesim']
#     lsAlgor = ['icebergAutosim', 'earlyStopAutosim', 'icebergEarlyStopAutosim']
#     lsAlgor = ['earlyStopAutosim']
#     lsAlgor = ['icebergEarlyStopAutosim']
#     lsAlgor = ['rolesim']
#     lsAlgor = ['simrank']
#     lsAlgor = ['autosim', 'vertBalAutosim']
    
        
    lllPartitions, numVerts = loadParts(lsPartFile, sMode)
    # convert to membership format
#     vPartMembership = [0 for x in range(numVerts)]
#     for (p, lPart) in enumerate(llPartitions):
# #         print len(lPart)
#         for v in lPart:
#             vPartMembership[v] = p

    for sAlgorname in lsAlgor:       
        print sAlgorname              
        fResultOut = open(sResultsOutPrefix + '_' + sAlgorname + '.out.csv' , 'a')
        csvResultOut = csv.writer(fResultOut, delimiter=',')               
                
        
        
        # store results
        lResults = []

        sTempOutFile = 'temp.out'
        sExec = '~/Programming/workspace/MuNS/cImpl/Release/vertSim'
        
        for dampingFactor in lLambdaRange:
            print [sExec, 
                                ' -t ' + str(maxIter), 
                                ' -e ' + str(convEpsilon), 
                                ' -i ' + sInitAlgor,
                                ' -b ' + str(ioBalance),
                                ' -d ' + str(dampingFactor),
                                ' -s ' + str(earlySimStopThres),
                                ' -c ' + str(icebergThres),
                                ' -a ' + str(iceBergAproxFactor),
                                sAlgorname, 
                                sTempOutFile,
                                ] + lsGraphFile
            proc = sp.Popen(sExec + 
                                ' -t ' + str(maxIter) + 
                                ' -e ' + str(convEpsilon) + 
                                ' -i ' + sInitAlgor +
                                ' -b ' + str(ioBalance) +
                                ' -d ' + str(dampingFactor) +
                                ' -s ' + str(earlySimStopThres) +
                                ' -c ' + str(icebergThres) +
                                ' -a ' + str(iceBergAproxFactor) +                                
                                ' ' + sAlgorname +
                                ' ' + sTempOutFile + ' ' +
                                ' '.join(lsGraphFile),
                                shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
            sOut, sErr = proc.communicate()
            print sOut
            print sErr
                
            # get the results
            with open(sTempOutFile, 'r') as fTempOutput:
                # read in first two rows
                sLine = fTempOutput.readline();
                lTokens = sLine.split(':')
                runningTime = float(lTokens[1])
                    
                sLine = fTempOutput.readline();
                lTokens = sLine.split(':')
                iterNum = int(lTokens[1])                        
                    
                # read in similarity matrix
                llSim = [];
                csvTempOut = csv.reader(fTempOutput, delimiter=',')
                for lRow in csvTempOut:
                    llSim.append([float(x) for x in lRow])
                
                lSim = [item for lSublist in llSim for item in lSublist]
                rowLen = len(llSim[0])                
                                            
                # do accuracy testing
                for (i, llPartitions1) in enumerate(lllPartitions):
                    assert(len(llPartitions1) > 0)

                    intra = avgIntraSim(llSim, llPartitions1)
                    inter = avgInterSim(llSim, llPartitions1)
                    intraInterRatio = intra / inter
                
                    intraRank, interRank = avgRank(lSim, rowLen, llPartitions1, simEpsilon)
                    
                    # average precision
    #                 precision, recall, fscore = avgRetrievalStats(llSim, llPartitions, topNCompared)
                    precisionAtK = avgPrecisionAtK(llSim, llPartitions1, topNCompared)
                    rPrecision = avgRPrecision(llSim, llPartitions1)
                    MAP = avgMAP(llSim, llPartitions1)
                    
                    currRow = [sAlgorname, sInitAlgor, ioBalance, dampingFactor, runningTime, iterNum, intra, inter, intraInterRatio, intraRank, interRank, precisionAtK, rPrecision, MAP]
                

                
                # store results
#                 lResults.append([sAlgorname, sInitAlgor, ioBalance, dampingFactor, runningTime, iterNum, intra, inter, intraInterRatio, intraRank, interRank, precisionAtK, rPrecision, MAP])
                
                    csvResultOut.writerow(currRow)
                    fResultOut.flush()
                    
                    for j in range(i+1, len(lllPartitions)):
                        llPartitions2 = lllPartitions[j];
                        
                        crossInter = avgCrossSim(llSim, llPartitions1, llPartitions2)
                        
                        crossIntraRank, crossInterRank = avgCrossRank(lSim, rowLen, llPartitions1, llPartitions2, simEpsilon)
                        
                        currRow = [i, i+1, crossInter, crossIntraRank, crossInterRank]
                        csvResultOut.writerow(currRow)
                        fResultOut.flush()
#                         crossRank(lSim, rowLen, llPartition1, llPartition2, simEpsilon)

                                    
            
            # delete temp.out
            os.remove('temp.out')
                  
        fResultOut.close()
    
        

######################################################################

def usage(sProgname):
    print >> sys.stderr, sProgname + "[original gml file that is to be aligned] [gml file where we want to align to] [output of aligned gml file]"
    sys.exit(2)


if __name__ == '__main__':
    main()