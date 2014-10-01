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
        elif opt == '-i':
            bIceberg
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
#     lsAlgor = ['icebergAutosim']   
    lsAlgor = ['icebergEarlyStopAutosim']
    defaultIcebergThes = 0.01
    defaultAlpha = 0.5
    defaultEarlyStopThres = 0.01
    
    lEarlyStopThres = [0.005, 0.01, 0.03, 0.05, 0.1]
    lIcebergThres = [0.95, 0.9, 0.8, 0.7]
    lIcebergAlpha = [0.1, 0.3, 0.5, 0.7, 0.9]
    
    
    lSizes = [500, 1000, 2000, 5000]
#     lSizes = [1000, 2000, 5000]
#     lSizes = [500]



#     sTempOutFile = 'temp.out'
    sExec = '~/Programming/workspace/cVertSim/cImpl/Release/cVertsim'
    
    

    for graphSize in lSizes:
        sCurrGraphDir = os.path.join(sGraphDir, str(graphSize))
        lsFiles = [ f for f in os.listdir(sCurrGraphDir) if os.path.isfile(os.path.join(sCurrGraphDir,f)) ]
                
        # store output
        hllResults = {}
        for sAlgorname in lsAlgor:
            # initialise hllResults according to the algorithm tested
            if sAlgorname == 'earlyStopAutosim':
                hllResults[sAlgorname] = [[] for x in range(len(lEarlyStopThres))]
                for p in range(len(lEarlyStopThres)):
                     hllResults[sAlgorname][p] = [[] for y in lsFiles]            
                
            elif sAlgorname == 'icebergAutosim':
                hllResults[sAlgorname] = [[] for x in range(len(lIcebergThres))]
                for t in range(len(lIcebergThres)):
                     hllResults[sAlgorname][t] = [[] for y in range(len(lIcebergAlpha))]
                     for a in range(len(lIcebergAlpha)):
                         hllResults[sAlgorname][t][a] = [[] for y in lsFiles]
                
            elif sAlgorname == 'icebergEarlyStopAutosim':
                hllResults[sAlgorname] = [[] for x in range(len(lEarlyStopThres))]                        
                for p in range(len(lEarlyStopThres)):
                    hllResults[sAlgorname][p] = [[] for y in range(len(lIcebergThres))]   
                    for t in range(len(lIcebergThres)):
                         hllResults[sAlgorname][p][t] = [[] for y in range(len(lIcebergAlpha))]
                         for a in range(len(lIcebergAlpha)):
                             hllResults[sAlgorname][p][t][a] = [[] for y in lsFiles]            

    
        
        # loop through the graph files
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
                sLine = fOrig.readline()
                lTokens = sLine.split(':')
                runningTime = float(lTokens[1])
                    
                sLine = fOrig.readline()
                lTokens = sLine.split(':')
                origIterNum = int(lTokens[1])                        
                    
                # read in similarity matrix
                csvOrig = csv.reader(fOrig, delimiter=',')
                for lRow in csvOrig:
                    lOrigSim.extend([float(x) for x in lRow])
                            
            assert(len(lOrigSim) > 0)
            
            
            for (i,sAlgorname) in enumerate(lsAlgor):
                print sAlgorname
                
                if sAlgorname == 'earlyStopAutosim':
                    testParas(sExec, hllResults, g, sAlgorname,  origIterNum, sInitAlgor, dampingFactor, sCurrGraphDir, sGraphFile[0:-extLen], sGraphFileExt, sResultsOutPrefix, defaultIcebergThes, defaultAlpha, defaultEarlyStopThres, lOrigSim, lEarlyStopThres)
                elif sAlgorname == 'icebergAutosim':
                    testParas(sExec, hllResults, g, sAlgorname,  origIterNum, sInitAlgor, dampingFactor, sCurrGraphDir, sGraphFile[0:-extLen], sGraphFileExt, sResultsOutPrefix, defaultIcebergThes, defaultAlpha, defaultEarlyStopThres, lOrigSim, lIcebergThres, lIcebergAlpha)
                elif sAlgorname == 'icebergEarlyStopAutosim':
                    testParas(sExec, hllResults, g, sAlgorname,  origIterNum, sInitAlgor, dampingFactor, sCurrGraphDir, sGraphFile[0:-extLen], sGraphFileExt, sResultsOutPrefix, defaultIcebergThes, defaultAlpha, defaultEarlyStopThres, lOrigSim, lEarlyStopThres, lIcebergThres, lIcebergAlpha)
 
        # write out results
        with open(sResultsOutPrefix + '_' + str(graphSize) + '.out.csv' , 'w') as fResultOut:
            csvResultOut = csv.writer(fResultOut, delimiter=',')
            for (sAlgor, lParaRunResults) in hllResults.items():
                
                if sAlgorname == 'earlyStopAutosim':
                    for (p, lRunResult) in enumerate(lParaRunResults):
                        lOut = [sAlgor, lEarlyStopThres[p]]
                        for (g, lResult) in enumerate(lRunResult):
                            lOut2 = lOut + [g] + lResult      
                            csvResultOut.writerow(lOut2)         
                        
                elif sAlgorname == 'icebergAutosim':
                    for (t, lThresResults) in enumerate(lParaRunResults):
                         for (a,lRunResult) in enumerate(lThresResults):
                             lOut = [sAlgor, lIcebergThres[t], lIcebergAlpha[a]]
                             for (g, lResult) in enumerate(lRunResult):
                                 lOut2 = lOut + [g] + lResult      
                                 csvResultOut.writerow(lOut2)    
                             
                elif sAlgorname == 'icebergEarlyStopAutosim':      
                    for (p, lESResult) in enumerate(lParaRunResults):
                        for (t, lThresResults) in enumerate(lESResult):
                             for (a,lRunResult) in enumerate(lThresResults):
                                 lOut = [sAlgor, lEarlyStopThres[p], lIcebergThres[t], lIcebergAlpha[a]]
                                 for (g, lResult) in enumerate(lRunResult):
                                     lOut2 = lOut + [g] + lResult      
                                     csvResultOut.writerow(lOut2)   
                
            
    #             lOut = [sAlgor]
                


######################################################################

def testParas(sExec, hllResults, graphIndex, sAlgorname, origIterNum, sInitAlgor, dampingFactor, sCurrGraphDir, sGraphFile, sExt, sResultsOutPrefix,   defaultIcebergThes, defaultAlpha, defaultEarlyStopThres, lOrigSim, *args):
    

    
    if sAlgorname == 'earlyStopAutosim':
        assert(len(args) == 1)
        lThres = args[0]
        
        for (p,para) in enumerate(lThres):
            sTempOutFile = sResultsOutPrefix + '_' + sGraphFile + '_' + sAlgorname + '_' + sInitAlgor + '_' + str(defaultIcebergThes) + '_' + str(defaultAlpha) + '_' + str(para) +  '.result.csv'
                                                                          
            runAlgor(sExec, origIterNum, sInitAlgor, dampingFactor, [('-s', para)], sCurrGraphDir, sGraphFile + sExt, sAlgorname, sTempOutFile)
            
            runningTime, iterNum, totalDiff = compare(sTempOutFile, lOrigSim)
    
            # store results
#             llResults.append([runningTime, iterNum, totalDiff])
            
            hllResults[sAlgorname][p][graphIndex] = [runningTime, iterNum, totalDiff]
                            
            os.remove(sTempOutFile)            
        
    elif sAlgorname == 'icebergAutosim':
        print args
        assert(len(args) == 2)
        lThres = args[0]
        lAlpha = args[1]
        
        for (t,thres) in enumerate(lThres):
            for (a, alpha) in enumerate(lAlpha):
                sTempOutFile = sResultsOutPrefix + '_' + sGraphFile + '_' + sAlgorname + '_' + sInitAlgor + '_' + str(thres) + '_' + str(alpha) + '_' + str(defaultEarlyStopThres) +  '.result.csv'
                                                                              
                runAlgor(sExec, origIterNum, sInitAlgor, dampingFactor, [('-c', thres), ('-a', alpha)], sCurrGraphDir, sGraphFile + sExt, sAlgorname, sTempOutFile)
                
                runningTime, iterNum, totalDiff = compare(sTempOutFile, lOrigSim)
        
                # store results
    #             llResults.append([runningTime, iterNum, totalDiff])
                
                hllResults[sAlgorname][t][a][graphIndex] = [runningTime, iterNum, totalDiff]
                                
                os.remove(sTempOutFile)           
        
        
    elif sAlgorname == 'icebergEarlyStopAutosim':
        assert(len(args) == 3)
        lESThres = args[0]
        lIBThres = args[1]
        lIBAlpha = args[2]
        
        for (p,esThres) in enumerate(lESThres):
            for (t,ibThres) in enumerate(lIBThres):
                for (a, ibAlpha) in enumerate(lIBAlpha):
                    sTempOutFile = sResultsOutPrefix + '_' + sGraphFile + '_' + sAlgorname + '_' + sInitAlgor + '_' + str(ibThres) + '_' + str(ibAlpha) + '_' + str(esThres) +  '.result.csv'
                                                                                  
                    runAlgor(sExec, origIterNum, sInitAlgor, dampingFactor, [('-s', esThres), ('-c', ibThres), ('-a', ibAlpha)], sCurrGraphDir, sGraphFile + sExt, sAlgorname, sTempOutFile)
                    
                    runningTime, iterNum, totalDiff = compare(sTempOutFile, lOrigSim)
            
                    # store results
        #             llResults.append([runningTime, iterNum, totalDiff])
                    
                    hllResults[sAlgorname][p][t][a][graphIndex] = [runningTime, iterNum, totalDiff]
                                    
                    os.remove(sTempOutFile)          
        
    else:
        pass
    
    




def runAlgor(sExec, origIterNum, sInitAlgor, dampingFactor, ltPara, sCurrGraphDir, sGraphFile, sAlgorname, sTempOutFile):
    """
    Run the algorihm.
    """
    
    print ltPara
    
    print [sExec, 
                    ' -t ' + str(origIterNum), 
                    ' -i ' + sInitAlgor,
                    ' -d ' + str(dampingFactor)] +\
                    [' {0} {1} '.format(option, para) for (option, para) in ltPara] +\
                    [
                    os.path.join(sCurrGraphDir, sGraphFile),
                    sAlgorname, 
                    sTempOutFile
                    ]
                    
    proc = sp.Popen(sExec + 
                        ' -t ' + str(origIterNum) + 
                        ' -i ' + sInitAlgor +
    #                                 ' -b ' + str(ioBalance) +
                        ' -d ' + str(dampingFactor) +
                        ' '.join([' {0} {1} '.format(option, para) for (option, para) in ltPara]) +
                        ' ' + os.path.join(sCurrGraphDir, sGraphFile) +
                        ' ' + sAlgorname +
                        ' ' + sTempOutFile
                        ,
                        shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
    
    sOut, sErr = proc.communicate()
    print sOut
    print sErr   

    


def compare(sTempOutFile, lOrigSim):

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
            
        return runningTime, iterNum, totalDiff    



def usage(sProgname):
    print >> sys.stderr, sProgname + "[original gml file that is to be aligned] [gml file where we want to align to] [output of aligned gml file]"
    sys.exit(2)


if __name__ == '__main__':
    main()