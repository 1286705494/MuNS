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
    sResultsOutPrefix = remainArgs[2]
    bApprox = int(remainArgs[3])

    
    
    
#     lLambdaRange = [0.1, 0.3, 0.5, 0.7, 0.9]
    
#     lBetaRange = [0.5]
#     lLambdaRange = [0.9]
    sInitAlgor = 'degBinaryInit'
#     lsInit = ['degRatioInit']
#     sOrigAlgor = 'autosim'
    lsAlgor = ['autosim']   
#     lsAlgor = ['icebergAutosim', 'earlyStopAutosim', 'icebergEarlyStopAutosim']
    defaultIcebergThes = 0.01
    defaultAlpha = 0.5
    defaultEarlyStopThres = 0.01
    llApproxPara = [[defaultIcebergThes, defaultAlpha, defaultEarlyStopThres], [defaultIcebergThes, defaultAlpha, defaultEarlyStopThres], [defaultIcebergThes, defaultAlpha, defaultEarlyStopThres]]

# 
#     sTempOutFile = 'temp.out'
    sExec = '~/Programming/workspace/cVertSim/cImpl/Release/cVertsim'
    
    lsFiles = [ f for f in os.listdir(sGraphDir) if os.path.isfile(os.path.join(sGraphDir,f)) ]
    
    extLen =  len(sGraphFileExt)
    for sGraphFile in lsFiles:
        if sGraphFile[-extLen:] != sGraphFileExt:
            continue
        
        for (i,sAlgorname) in enumerate(lsAlgor):       
            print sAlgorname           
            print sGraphFile[0:-extLen]
            if bApprox > 0:   
                lPara = llApproxPara[i]
                sTempOutFile = sResultsOutPrefix + '_' + sGraphFile[0:-extLen] + '_' + sAlgorname + '_' + sInitAlgor + '_' + str(lPara[0]) + '_' + str(lPara[1]) + '_' + str(lPara[2]) +  '.result.csv'
            else:
                sTempOutFile = sResultsOutPrefix + '_' + sGraphFile[0:-extLen] + '_' + sAlgorname + '_' + sInitAlgor + '.result.csv'                                              
                
            # store results
            lResults = []
    
            print [sExec, 
                                ' -t ' + str(maxIter), 
                                ' -e ' + str(convEpsilon), 
                                ' -i ' + sInitAlgor,
#                                 ' -b ' + str(ioBalance),
                                ' -d ' + str(dampingFactor),
                                sGraphFile,
                                sAlgorname, 
                                sTempOutFile
                                ]
            proc = sp.Popen(sExec + 
                                ' -t ' + str(maxIter) + 
                                ' -e ' + str(convEpsilon) + 
                                ' -i ' + sInitAlgor +
#                                 ' -b ' + str(ioBalance) +
                                ' -d ' + str(dampingFactor) +
                                ' ' + os.path.join(sGraphDir, sGraphFile) +
                                ' ' + sAlgorname +
                                ' ' + sTempOutFile
                                ,
                                shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
            sOut, sErr = proc.communicate()
            print sOut
            print sErr        
    

                

    
        

######################################################################

def usage(sProgname):
    print >> sys.stderr, sProgname + "[original gml file that is to be aligned] [gml file where we want to align to] [output of aligned gml file]"
    sys.exit(2)


if __name__ == '__main__':
    main()