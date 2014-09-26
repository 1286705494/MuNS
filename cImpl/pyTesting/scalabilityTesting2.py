'''
Created on 26/09/2014

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
        options = "t:e:d:"
        # get options
        optList, remainArgs = getopt.gnu_getopt(sys.argv[1:], options)
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
        usage(sys.argv[0])
    
      
    maxIter = 20
    convEpsilon = 0.01
    dampingFactor = 0.9
    for opt, arg in optList:
        if opt == '-t':
            maxIter = int(arg)     
        elif  opt == '-e':
            convEpsilon = float(arg)
        elif opt == '-d':
            dampingFactor = float(arg)
        else:
            print >> sys.stderr, sys.argv[0] + " -m <mat index wanted> [edge list]"
            sys.exit(2)    
    
    # check number of arguments
    if len(remainArgs) != 4:
        usage(sys.argv[0])
    
    sGraphDir = remainArgs[0]    
    sGraphFileExt = remainArgs[1]
    sResultsDir = remainArgs[2]
    sResultsOutPrefix = remainArgs[3]

    
    
#     lLambdaRange = [0.1, 0.3, 0.5, 0.7, 0.9]
    
#     lBetaRange = [0.5]
#     lLambdaRange = [0.9]
    sInitAlgor = 'degBinaryInit'
#     lsInit = ['degRatioInit']
    sOrigAlgor = 'autosim'   
    lsAlgor = ['earlyStopAutosim']
    defaultIcebergThes = 0.01
    defaultAlpha = 0.5
    defaultEarlyStopThres = 0.01
    
    lPara = [0.005, 0.01, 0.03, 0.05, 0.1]
    
    lSizes = [1000, 2000, 5000]



#     sTempOutFile = 'temp.out'
    sExec = '~/Programming/workspace/cVertSim/cImpl/Release/cVertsim'
    
    

    for graphSize in lSizes:
        sCurrGraphDir = os.path.join(sGraphDir, str(graphSize))
        lsFiles = [ f for f in os.listdir(sCurrGraphDir) if os.path.isfile(os.path.join(sCurrGraphDir,f)) ]
                
        # store output
        hllResults = {}
        for sAlgor in lsAlgor:
            hllResults[sAlgor] = [[] for x in lPara]
            
            for p in range(len(lPara)):
                 hllResults[sAlgor][p] = [[] for y in lsFiles]
    
        
        extLen =  len(sGraphFileExt)
        for (g,sGraphFile) in enumerate(lsFiles):
            if sGraphFile[-extLen:] != sGraphFileExt:
                continue
            
            print sGraphFile
            
            # load the original results
            sCurrResultsDir = os.path.join(sResultsDir, str(graphSize))
            lResultsFiles = glob.glob(os.path.join(sCurrResultsDir, '*_' + sGraphFile[0:-extLen] + '_*'))
            assert(len(lResultsFiles) == 1)
            sOrigResultFile = lResultsFiles[0]
            
            print sOrigResultFile
            
            lOrigSim = [];
            origIterNum = 20
            with open(sOrigResultFile, 'r') as fOrig:
                # read in first two rows
                sLine = fOrig.readline();
                lTokens = sLine.split(':')
                runningTime = float(lTokens[1])
                    
                sLine = fOrig.readline();
                lTokens = sLine.split(':')
                origIterNum = int(lTokens[1])                        
                    
                # read in similarity matrix
                
                csvOrig = csv.reader(fOrig, delimiter=',')
                for lRow in csvOrig:
                    lOrigSim.extend([float(x) for x in lRow])
                            
            assert(len(lOrigSim) > 0)
            
            
            for (i,sAlgorname) in enumerate(lsAlgor):       
                print sAlgorname           
                
                for (p,para) in enumerate(lPara):
                    sTempOutFile = sResultsOutPrefix + '_' + sGraphFile[0:-extLen] + '_' + sAlgorname + '_' + sInitAlgor + '_' + str(defaultIcebergThes) + '_' + str(defaultAlpha) + '_' + str(para) +  '.result.csv'
                                                                              
                    # store results
                    lResults = []
            
                    print [sExec, 
                                        ' -t ' + str(origIterNum), 
                                        ' -i ' + sInitAlgor,
        #                                 ' -b ' + str(ioBalance),
                                        ' -d ' + str(dampingFactor),
                                        ' -s ' + str(para),
                                        sGraphFile,
                                        sAlgorname, 
                                        sTempOutFile
                                        ]
                    proc = sp.Popen(sExec + 
                                        ' -t ' + str(origIterNum) + 
                                        ' -i ' + sInitAlgor +
        #                                 ' -b ' + str(ioBalance) +
                                        ' -d ' + str(dampingFactor) +
                                        ' -s ' + str(para) +
                                        ' ' + os.path.join(sCurrGraphDir, sGraphFile) +
                                        ' ' + sAlgorname +
                                        ' ' + sTempOutFile
                                        ,
                                        shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
                    sOut, sErr = proc.communicate()
                    print sOut
                    print sErr        
        
                    # load file then compare with existing
    
                    with open(sTempOutFile, 'r') as fTempOutput:
                        # read in first two rows
                        sLine = fTempOutput.readline();
                        lTokens = sLine.split(':')
                        runningTime = float(lTokens[1])
                            
                        sLine = fTempOutput.readline();
                        lTokens = sLine.split(':')
                        iterNum = int(lTokens[1])                        
                            
                        # read in similarity matrix
                        lCurrSim = [];
                        csvTempOut = csv.reader(fTempOutput, delimiter=',')
                        for lRow in csvTempOut:
                            lCurrSim.extend([float(x) for x in lRow])                
    
                        assert(len(lOrigSim) == len(lCurrSim))
                        # compare the outputs
                        totalDiff = 0;
                        for i in range(len(lOrigSim)):
                            totalDiff += abs(lOrigSim[i] - lCurrSim[i])
                            
                        # store results
                        hllResults[sAlgorname][p][g] = [runningTime, iterNum, totalDiff]
                                    
                    
                    os.remove(sTempOutFile)
        
        # write out results
        with open(sResultsOutPrefix + '_' + str(graphSize) + '_' + '.out.csv' , 'w') as fResultOut:
            csvResultOut = csv.writer(fResultOut, delimiter=',')
            for (sAlgor, lParaRunResults) in hllResults.items():
    #             lOut = [sAlgor]
                for (p, lRunResult) in enumerate(lParaRunResults):
                    lOut = [sAlgor, lPara[p]]
                    for (g, lResult) in enumerate(lRunResult):
                        lOut2 = lOut + [g] + lResult      
                        csvResultOut.writerow(lOut2)

######################################################################

def usage(sProgname):
    print >> sys.stderr, sProgname + "[original gml file that is to be aligned] [gml file where we want to align to] [output of aligned gml file]"
    sys.exit(2)


if __name__ == '__main__':
    main()