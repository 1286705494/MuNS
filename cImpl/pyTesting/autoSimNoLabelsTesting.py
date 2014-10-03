'''
Created on 10/09/2014

@author: Jeffrey Chan, 2014
'''
    
import sys
import getopt
import glob
import csv
import os.path
import string
import subprocess as sp
import numpy as np
from testingUtils import *
from kmedriods import cluster
import networkx as nx

def main():
    # process options
    try:
        # option list
        options = "t:c:b:k:m:se:o"
        # get options
        optList, remainArgs = getopt.gnu_getopt(sys.argv[1:], options)
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
        usage(sys.argv[0])
    
     
    maxIter = 20
    convEpsilon = 0.03
    ioBalance = 0.5
    clusNum = 2
    inBinNum = 10
    outBinNum = 10
    bShellAnalysis = False
    simEpsilon = 0.001
    bVertSubtractOne = False
    for opt, arg in optList:
        if opt == '-t':
            maxIter = int(arg)
        elif  opt == '-c':
            convEpsilon = float(arg)
        elif  opt == '-b':
            ioBalance = float(arg)
        elif  opt == '-k':
            clusNum = int(arg)
        elif  opt == '-m':
            lBin = string.split(arg, ',')
            assert(len(lBin) == 2)
            outBinNum = int(lBin[0])
            inBinNum = int(lBin[1])            
        elif opt == "-s":
            bShellAnalysis = True
        elif  opt == '-e':
            simEpsilon = float(arg)
            assert(simEpsilon <= 1.0 and simEpsilon >= 0.0)      
        elif opt == '-o':
            bVertSubtractOne = True         
        else:
            print >> sys.stderr, sys.argv[0] + " -m <mat index wanted> [edge list]"
            sys.exit(2)    
    
    # check number of arguments
    if len(remainArgs) != 2:
        usage(sys.argv[0])
    
    sGraphFile = remainArgs[0]    
    sResultsOutPrefix = remainArgs[1]
    
    
    
#     lLambdaRange = [0.1, 0.3, 0.5, 0.7, 0.9]
    lLambdaRange = [0.9]
#     lBetaRange = [0.5]
#     lLambdaRange = [0.9]
#     sInitAlgor = 'degBinaryInit'
    sInitAlgor = 'binaryInit'
#     lsInit = ['degRatioInit']
    lsAlgor = ['autosim', 'simrank', 'prank', 'rolesim', 'matchsim']
#     lsAlgor = ['simrank']
#     lsAlgor = ['autosim']
    
    

    for sAlgorname in lsAlgor:       
        print sAlgorname              
        fResultOut = open(sResultsOutPrefix + '_' + sAlgorname + '.out.csv' , 'w')
        csvResultOut = csv.writer(fResultOut, delimiter=',')               
                        
        # store results
        lResults = []

        sTempOutFile = 'temp.out'
#         sExec = '~/Programming/MuNS/cImpl/Release/cVertSim'
        sExec = '~/Programming/workspace/cVertSim/cImpl/Release/cVertsim'
        
        for dampingFactor in lLambdaRange:
            lPara = [
                                ' -t ' + str(maxIter), 
                                ' -e ' + str(convEpsilon), 
                                ' -i ' + sInitAlgor,
#                                 ' -b ' + str(ioBalance),
                                ' -d ' + str(dampingFactor)
                     ]
            
            if bVertSubtractOne:
                lPara.append(' -m ')
                
            lExecString = [sExec] + lPara + [sGraphFile,sAlgorname, sTempOutFile] 
            print lExecString
                                
            proc = sp.Popen(" ".join(lExecString),
                                shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
            sOut, sErr = proc.communicate()
            print sOut
            print sErr
                
            # get the results
            with open('temp.out', 'r') as fTempOutput:
                # read in first two rows
                sLine = fTempOutput.readline();
                lTokens = sLine.split(':')
                runningTime = float(lTokens[1])
                    
                sLine = fTempOutput.readline();
                lTokens = sLine.split(':')
                iterNum = int(lTokens[1])                        
                    
                # read in similarity matrix
                llDis = [];
                llSim = [];
                csvTempOut = csv.reader(fTempOutput, delimiter=',')
                for lRow in csvTempOut:
                    llDis.append([1-float(x) for x in lRow])
                    if bShellAnalysis:
                        llSim.append([float(x) for x in lRow])
                                            
                # perform k-mediods clustering
                mDis = np.array(llDis)
                vClusters, vMedriods = cluster(mDis, clusNum)
                
                # compute compactness
                compactNess = computeCompactness(mDis, vClusters, vMedriods)
                
                # store results
                lResults.append([sAlgorname, sInitAlgor, ioBalance, dampingFactor, runningTime, iterNum, compactNess])
                currRow = [sAlgorname, sInitAlgor, ioBalance, dampingFactor, runningTime, iterNum, compactNess]                
                
                if bShellAnalysis:
                    # compute d-shell
                    with open(sGraphFile, 'r') as fGraph:
                        csvGraph = csv.reader(fGraph, delimiter=',')
                        # read into network graph
                        gGraph = nx.DiGraph()
                        for lRow in csvGraph:
#                             print lRow
                            src = int(lRow[0])
                            tar = int(lRow[1])
                            
                            if bVertSubtractOne:
                                src -= 1
                                tar -= 1
                            
                            gGraph.add_edge(src, tar)
                        
                        llDShells = computeDShell(gGraph, outBinNum, inBinNum)
                        
                        lSim = [item for lSublist in llSim for item in lSublist]
                        rowLen = len(llSim[0])
#                         intraRank, interRank = avgRank(lSim, rowLen, llDShells, simEpsilon)
                        
                        lRank = percentileRank(lSim, simEpsilon)
                        mIntraRank, mShellInterRank = intraCrossShellsRank(llDShells, lRank, outBinNum, inBinNum)
                        
                        print mIntraRank
                        print mShellInterRank
                        
#                     lSim = [item for lSublist in llSim for item in lSublist]
#                     lHist = np.histogram(lSim, bins=10)
                    
                        
#                     lResults[-1].extend([intraRank, interRank])
#                     currRow.extend([intraRank, interRank])

                csvResultOut.writerow(currRow)
                fResultOut.flush()
                
            
            # delete temp.out
            os.remove('temp.out')
                  
        fResultOut.close()
    
        

######################################################################

def usage(sProgname):
    print >> sys.stderr, sProgname + "[original gml file that is to be aligned] [gml file where we want to align to] [output of aligned gml file]"
    sys.exit(2)


if __name__ == '__main__':
    main()