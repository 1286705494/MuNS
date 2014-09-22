'''
Created on 07/09/2014

@author: Jeffrey Chan, 2014
'''

import sys
import getopt
import csv
import os.path
import string
import networkx as nx


def main():
    
    # process options
    try:
        # option list
        options = "c"
        # get options
        optList, remainArgs = getopt.gnu_getopt(sys.argv[1:], options)
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
        usage(sys.argv[0])
        
        
    bConnComp = False
    for opt, arg in optList:
        if opt == "-c":
            bConnComp = True
        else:
            print >> sys.stderr, sys.argv[0] + " -m <mat index wanted> [edge list]"
            sys.exit(2)    
    
        # check number of arguments
    if len(remainArgs) != 3:
        usage(sys.argv[0])
    
    sCiteFile = remainArgs[0]    
    sMembershipFile = remainArgs[1]
    sOutPrefix = remainArgs[2]
    
    hMembership = {}
    hUniqClass = {}
    hVertMap = {}
    currMembClassNum = 0
    currVertId = 0    

    gCite = nx.DiGraph()

    with open(sMembershipFile, 'r') as fMembership:
        csvMembership = csv.reader(fMembership, delimiter=',')
        for lsRow in csvMembership:
            vertId = lsRow[0]
            sMembershipClass = lsRow[-1]
            hVertMap[vertId] = currVertId;
            
            
            gCite.add_node(vertId, membership=sMembershipClass)
            
            hMembership[vertId] = sMembershipClass
            if sMembershipClass not in hUniqClass:
                hUniqClass[sMembershipClass] = currMembClassNum;
                currMembClassNum += 1
            
            currVertId += 1      
  

    with open(sCiteFile, 'r') as fCite:
        for sLine in fCite:
            lTokens = string.split(string.strip(sLine), ',')
            src = lTokens[0]
            tar = lTokens[1]
            if src in hVertMap and tar in hVertMap:
                gCite.add_edge(lTokens[0], lTokens[1])
                
        if bConnComp:
            llConnComp = sorted(nx.weakly_connected_components(gCite), key = len, reverse=True)

            for comp in range(1, len(llConnComp)):
                gCite.remove_nodes_from(llConnComp[comp])
                for v in llConnComp[comp]:
                    del hVertMap[v]
                    del hMembership[v]
                
            # produce new mapping 
            hNewVertMap = {}
            currVertId = 0
            for origId in hVertMap.keys():
                hNewVertMap[origId] = currVertId
                currVertId += 1
            hVertMap = hNewVertMap
                
                
    
        print gCite.number_of_nodes()
        print len(gCite.edges())
        

    # write out as vertId and integer membership
    with open(sOutPrefix + 'Membership.csv', 'w') as fMembershipOut:
        csvMembOut = csv.writer(fMembershipOut, delimiter=',')
        for (vertId, membershipClass) in hMembership.items():
            csvMembOut.writerow([hVertMap[vertId], hUniqClass[membershipClass]])
            
    # count number of vertices with classes and only write out edges that have a class attribute
    with open(sOutPrefix + 'Adj.csv', 'w') as fCiteOut:
        csvCiteOut = csv.writer(fCiteOut, delimiter=',')    
        for (src, tar) in gCite.edges_iter():
            if src in hMembership and tar in hMembership:
                csvCiteOut.writerow([hVertMap[src], hVertMap[tar]])
                
                
    with open(sOutPrefix + 'Map.csv', 'w') as fMapOut:
        csvMapOut = csv.writer(fMapOut, delimiter=',')    
        for (orig, mapped) in hVertMap.items():
            csvMapOut.writerow([mapped, orig]) 
            

            
########################################################################


def usage(sProgname):
    print >> sys.stderr, sProgname + "[original gml file that is to be aligned] [gml file where we want to align to] [output of aligned gml file]"
    sys.exit(2)


if __name__ == '__main__':
    main()