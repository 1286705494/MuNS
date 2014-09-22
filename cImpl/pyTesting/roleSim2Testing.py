'''
Created on 15/08/2014

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

def main():
    # process options
    try:
        # option list
        options = "p:"
        # get options
        optList, remainArgs = getopt.gnu_getopt(sys.argv[1:], options)
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
        usage(sys.argv[0])
    
     
    sMode = 'list'
    for opt, arg in optList:
        if opt == "-p":
            sMode = arg
        else:
            print >> sys.stderr, sys.argv[0] + " -m <mat index wanted> [edge list]"
            sys.exit(2)    
    
    # check number of arguments
    if len(remainArgs) != 5:
        usage(sys.argv[0])
    
    sBaseDir = remainArgs[0]    
    sFilePrefix = remainArgs[1]
    sFileGraphSuffix = remainArgs[2]
    sFilePartSuffix = remainArgs[3]
#     vertNum = remainArgs[4]
    sResultsOutPrefix = remainArgs[4]
    
    
    lBetaRange = [0, 0.25, 0.5, 0.75, 1]
    lLambdaRange = [0.1, 0.3, 0.5, 0.7, 0.9]
#     lBetaRange = [0.5]
#     lLambdaRange = [0.9]
    lsInit = ['binaryInit', 'degBinaryInit', 'degRatioInit']
#     lsInit = ['degRatioInit']
    sAlgor = 'autosim'
    
    # subdirs
    lsSubDirs = ['hierarchy', 'corePeriphery']
#     lsSizeDirs = ['500', '1000', '2000']
    lsSizeDirs = ['500']
    
    for sSubDir in lsSubDirs:
        for sSizeDir in lsSizeDirs:
    
            vertNum = int(sSizeDir)
            # get the graphs 
            lsFiles = glob.glob(os.path.join(sBaseDir, sSubDir, sSizeDir, '*' + sFilePrefix + '*' + sFileGraphSuffix));
            
            currFileId = 1
            
            fResultOut = open(sResultsOutPrefix + '_' + sSubDir + '_' + sSizeDir + '.out.csv' , 'w')
            csvResultOut = csv.writer(fResultOut, delimiter=',')       
            
            for sGraphFile in lsFiles:
                print sGraphFile
                # make sure we have corresponding part file as well
                index = string.rfind(sGraphFile, sFileGraphSuffix)
                sPartFile = sGraphFile[0:index] + sFilePartSuffix
                # make them absolute
                sAbsGraphFile = os.path.abspath(sGraphFile)
                sAbsPartFile = os.path.abspath(sPartFile)
        
                
                bFileExist = os.path.isfile(sAbsPartFile);
                if not bFileExist:
                    print >> sys.stderr, "Cannot locate partition file " + sAbsPartFile 
                    usage(sys.argv[0])
                     
                # run each graph
        #         print sAbsGraphFile
        #         print sAbsPartFile
                
#                 llPartitions = []
                
                # store results
                lResults = []
                
                sTempOutFile = 'temp.out'
                sAlgorname = 'autosim'
                maxIter = 20
                convEpsilon = 0.03
                
                llPartitions = loadPart(sAbsPartFile, sMode)
                
 
                        
        #             print llPartitions
                
                
                for sInitAlgor in lsInit:
                    for ioBalance in lBetaRange:
                        for conv in lLambdaRange:
                            sExec = '~/Programming/workspace/matlabSimrank/rolesim2/cRolesim2/Release/cRoleSim2' 
                            print [sExec, 
                                                ' -t ' + str(maxIter), 
                                                ' -e ' + str(convEpsilon), 
                                                ' -i ' + sInitAlgor,
                                                ' -b ' + str(ioBalance),
                                                ' -d ' + str(conv),
                                                sAbsGraphFile,
                                                sAlgorname, 
                                                str(vertNum),
                                                sTempOutFile
                                                ]
                            proc = sp.Popen(sExec + 
                                                ' -t ' + str(maxIter) + 
                                                ' -e ' + str(convEpsilon) + 
                                                ' -i ' + sInitAlgor +
                                                ' -b ' + str(ioBalance) +
                                                ' -d ' + str(conv) +
                                                ' ' + sAbsGraphFile +
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
                                    
        #                         print llSim
                                        
                                # do accuracy testing
                                assert(len(llPartitions) > 0)
                                intra = avgIntraSim(llSim, llPartitions)
                                inter = avgInterSim(llSim, llPartitions)
                                intraInterRatio = intra / inter
                                    
                                # store results
                                lResults.append([sInitAlgor, ioBalance, conv, runningTime, iterNum, intra, inter, intraInterRatio])
                                currRow = [currFileId]
                                currRow.extend([sInitAlgor, ioBalance, conv, runningTime, iterNum, intra, inter, intraInterRatio])
                                csvResultOut.writerow(currRow)
                                fResultOut.flush()
                                
                            
                            # delete temp.out
                            os.remove('temp.out')
                              
                             
#                 for lRow in lResults:
#                     currRow = [currFileId]
#                     currRow.extend(lRow)
#                     csvResultOut.writerow(currRow)

                
                    
                                
                                
                currFileId += 1
                    
            fResultOut.close()
         
        




######################################################################

def usage(sProgname):
    print >> sys.stderr, sProgname + "[original gml file that is to be aligned] [gml file where we want to align to] [output of aligned gml file]"
    sys.exit(2)


if __name__ == '__main__':
    main()