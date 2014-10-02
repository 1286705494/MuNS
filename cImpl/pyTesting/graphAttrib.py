
import getopt
import sys
import networkx as nx
import numpy as np


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
    
      
#     maxIter = 20
#     convEpsilon = 0.01
#     dampingFactor = 0.9
#     for opt, arg in optList:
#         if opt == '-t':
#             maxIter = int(arg)     
#         elif  opt == '-e':
#             convEpsilon = float(arg)
#         elif opt == '-d':
#             dampingFactor = float(arg)
#         elif opt == '-i':
#             bIceberg
#         else:
#             print >> sys.stderr, sys.argv[0] + " -m <mat index wanted> [edge list]"
#             sys.exit(2)    
    
    # check number of arguments
    if len(remainArgs) >= 2:
        usage(sys.argv[0])
    
    sResultsFilename = remainArgs[0]
    lsGraphFiles = remainArgs[1:]
    
    
    for sGraphFile in lsGraphFiles:
        with open(sGraphFile, 'r') as fGraph:
            
            gGraph = nx.DiGraph()
            
            csvGraph = csv.reader(fGraph, delimiter=',')
            # read into network graph
            for lRow in csvGraph:
    
                src = int(lRow[0])
                tar = int(lRow[1])
                
                gGraph.add_edge(src, tar)    
        
            binNum = 10
            
            # compute in/out degree distribution
            lOutDeg = list(gGraph.out_degree().values())
            lInDeg = list(gGraph.in_degree().values())
            lOutDegHist = np.histogram(lOutDeg, bins=binNum)
            lInDegHist = np.histogram(lInDeg, bins=binNum)
            
            
            # compute vert/edge ratio
            vertEdgeRatio = float(gGraph.number_of_nodes()) / gGraph.number_of_edges()
            
            # compute # of connected components
            ccNum = nx.number_connected_components(gGraph)
            strongCCNum = nx.number_strongly_connected_components(gGraph)
            
            # compute distribution of vert/edge ratio
            lVertInOutDeg = [float(gGraph.in_degree(v)) / gGraph.out_degree(v) for v in gGraph.nodes_iter()]
            lVertIODegHist = np.histogram(lVertInOutDeg, bins=binNum)
            
            # compute clustering coefficient
            avgClusCoeff = nx.average_clustering(gGraph)
            
            # degree assoritivity
            degAssort = nx.degree_assortativity_coefficient(gGraph)
            
            print vertEdgeRatio, ccNum, strongCCNum, avgClusCoeff, degAssort, lOutDegHist, lInDegHist, lVertIODegHist
    

#############################################################3
    
if __name__ == "__main__":
    main()