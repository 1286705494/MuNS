'''
Created on 25/09/2014

@author: Jeffrey Chan, 2014
'''

import sys
import getopt
import glob
import csv
import os
import os.path
import string
import subprocess as sp
import numpy as np
from testingUtils import *




def main():
    # process options
    try:
        # option list
        options = ""
        # get options
        optList, remainArgs = getopt.gnu_getopt(sys.argv[1:], options)
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
        usage(sys.argv[0])
    
#      
# 
#     for opt, arg in optList:
#         if opt == '-x':
#             pass
#         else:
#             print >> sys.stderr, sys.argv[0] + " -m <mat index wanted> [edge list]"
#             sys.exit(2)    
    
    # check number of arguments
    if len(remainArgs) != 3:
        usage(sys.argv[0])
    
    sGraphDir = remainArgs[0]    
    sGraphFileExt = remainArgs[1]
    sResultsOutPrefix = remainArgs[2]
    
    
    
#     lLambdaRange = [0.1, 0.3, 0.5, 0.7, 0.9]
    dampingFactor = 0.9
#     lBetaRange = [0.5]
#     lLambdaRange = [0.9]
    sInitAlgor = 'degBinaryInit'
#     lsInit = ['degRatioInit']
    sOrigAlgor = 'autosim'
    lsApproxAlgor = []


    sTempOutFile = 'temp.out'
    sExec = '~/Programming/workspace/matlabSimrank/rolesim2/cRolesim2/Release/cRoleSim2'
    

    lsFiles = [ f for f in os.listdir(mypath) if os.path.isfile(os.path.join(sGraphDir,f)) ]
    
    # number of vertices
    vertNum = X
    
    extLen = - len(sGraphFileExt)
    for sFile in lsFiles:
        if sFile[extLen:] != sGraphFileExt:
            continue
        
                
        print sOrigAlgor
        fResultOut = open(sResultsOutPrefix + '_' + sOrigAlgor + '.out.csv' , 'w')
        csvResultOut = csv.writer(fResultOut, delimiter=',')
                     
        # compute sOrigAlgor first
        print [sExec, 
                            ' -t ' + str(maxIter), 
                            ' -e ' + str(convEpsilon), 
                            ' -i ' + sInitAlgor,
#                                 ' -b ' + str(ioBalance),
                            ' -d ' + str(dampingFactor),
                            sGraphFile,
                            sOrigAlgor, 
                            str(vertNum),
                            sTempOutFile
                            ]
        proc = sp.Popen(sExec + 
                            ' -t ' + str(maxIter) + 
                            ' -e ' + str(convEpsilon) + 
                            ' -i ' + sInitAlgor +
#                                 ' -b ' + str(ioBalance) +
                            ' -d ' + str(dampingFactor) +
                            ' ' + sGraphFile +
                            ' ' + sOrigAlgor +
                            ' ' + str(vertNum) +
                            ' ' + sTempOutFile
                            ,
                            shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
        sOut, sErr = proc.communicate()
        print sOut
        print sErr        
        
    

    for sAlgorname in lsAlgor:       
        print sAlgorname              
        fResultOut = open(sResultsOutPrefix + '_' + sAlgorname + '.out.csv' , 'w')
        csvResultOut = csv.writer(fResultOut, delimiter=',')               
                
        llPartitions = loadPart(sPartFile, sMode)
        
        # store results
        lResults = []


        
        for dampingFactor in lLambdaRange:
            print [sExec, 
                                ' -t ' + str(maxIter), 
                                ' -e ' + str(convEpsilon), 
                                ' -i ' + sInitAlgor,
#                                 ' -b ' + str(ioBalance),
                                ' -d ' + str(dampingFactor),
                                sGraphFile,
                                sAlgorname, 
                                str(vertNum),
                                sTempOutFile
                                ]
            proc = sp.Popen(sExec + 
                                ' -t ' + str(maxIter) + 
                                ' -e ' + str(convEpsilon) + 
                                ' -i ' + sInitAlgor +
#                                 ' -b ' + str(ioBalance) +
                                ' -d ' + str(dampingFactor) +
                                ' ' + sGraphFile +
                                ' ' + sAlgorname +
                                ' ' + str(vertNum) +
                                ' ' + sTempOutFile
                                ,
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
                llSim = [];
                csvTempOut = csv.reader(fTempOutput, delimiter=',')
                for lRow in csvTempOut:
                    llSim.append([float(x) for x in lRow])
                                            
                # do accuracy testing
                assert(len(llPartitions) > 0)
                intra = avgIntraSim(llSim, llPartitions)
                inter = avgInterSim(llSim, llPartitions)
                intraInterRatio = intra / inter
                
                intraRank, interRank = avgRank(llSim, llPartitions, simEpsilon)
                    
                # average precision
#                 precision, recall, fscore = avgRetrievalStats(llSim, llPartitions, topNCompared)
                precisionAtK = avgPrecisionAtK(llSim, llPartitions, topNCompared)
                rPrecision = avgRPrecision(llSim, llPartitions)
                MAP = avgMAP(llSim, llPartitions)
                
                # store results
                lResults.append([sAlgorname, sInitAlgor, ioBalance, dampingFactor, runningTime, iterNum, intra, inter, intraInterRatio, intraRank, interRank, precisionAtK, rPrecision, MAP])
                currRow = [sAlgorname, sInitAlgor, ioBalance, dampingFactor, runningTime, iterNum, intra, inter, intraInterRatio, intraRank, interRank, precisionAtK, rPrecision, MAP]
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