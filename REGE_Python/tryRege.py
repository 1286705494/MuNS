__author__ = 'Youhan'

import networkx as nx
import sys

numIteration = 1000 #100000

fh = open(sys.argv[1], "rb")
diGraph = nx.read_edgelist(fh, delimiter=" ", create_using=nx.DiGraph(), nodetype=int, data=False)
fh.close()
numVertices = diGraph.number_of_nodes()
listVertices = diGraph.nodes()

print "there are", numVertices, "nodes in the graph"

rege = [[[1.0 for col in xrange(numVertices)] for row in xrange(numVertices)], [[0.0 for col in xrange(numVertices)] for row in xrange(numVertices)]]
denominator = list(diGraph.degree().values())

for i in xrange(1, numIteration + 1):
    print "current iteration:", "%8d/%d" % (i, numIteration), "\r",
    for v1 in listVertices:
        for v2 in listVertices:
            if v1 == v2:
                rege[i%2][v1][v2] = 1.0
                break
            rege[i%2][v1][v2] = 0.0
            for v3 in diGraph.successors(v1):
                maxv4 = -1
                for v4 in diGraph.successors(v2):
                    if maxv4 == -1 or rege[(i-1)%2][v3][v4] > rege[(i-1)%2][v3][maxv4]:
                        maxv4 = v4
                if maxv4 <> -1:
                    rege[i%2][v1][v2] += rege[(i-1)%2][v3][maxv4]
            for v3 in diGraph.successors(v2):
                maxv4 = -1
                for v4 in diGraph.successors(v1):
                    if maxv4 == -1 or rege[(i-1)%2][v3][v4] > rege[(i-1)%2][v3][maxv4]:
                        maxv4 = v4
                if maxv4 <> -1:
                    rege[i%2][v1][v2] += rege[(i-1)%2][v3][maxv4]
            for v3 in diGraph.predecessors(v1):
                maxv4 = -1
                for v4 in diGraph.predecessors(v2):
                    if maxv4 == -1 or rege[(i-1)%2][v3][v4] > rege[(i-1)%2][v3][maxv4]:
                        maxv4 = v4
                if maxv4 <> -1:
                    rege[i%2][v1][v2] += rege[(i-1)%2][v3][maxv4]
            for v3 in diGraph.predecessors(v2):
                maxv4 = -1
                for v4 in diGraph.predecessors(v1):
                    if maxv4 == -1 or rege[(i-1)%2][v3][v4] > rege[(i-1)%2][v3][maxv4]:
                        maxv4 = v4
                if maxv4 <> -1:
                    rege[i%2][v1][v2] += rege[(i-1)%2][v3][maxv4]
            rege[i%2][v1][v2] /= denominator[v1] + denominator[v2]
            rege[i%2][v2][v1] = rege[i%2][v1][v2]
print "finished"

if True:
    for v1 in listVertices:
        for v2 in listVertices:
            print "%5.3f" % rege[numIteration%2][v1][v2],
        print ''
else:
    exactCount = 0
    for v1 in listVertices:
        for v2 in listVertices:
            if v2 > v1:
                break
            if rege[numIteration%2][v1][v2] > 0.999999 and v1 <> v2:
                exactCount += 1
    print exactCount, "pairs of nodes are exactly equivalent under", numIteration, "iterations"
