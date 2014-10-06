/*
 * randGraph.cpp
 *
 *  Created on: 06/10/2014
 *      Author: youhan
 */

#include <iostream>
#include <fstream>


#include "randGraph.h"

using namespace std;

void unifDistGraph(string fileName, int vertNum, int edgeNum /* = 0 */, bool nonMultiEdge /* = true */)
{
	srand (time(NULL));

	ofstream fout(fileName);

	int* mAdjacent = new int[vertNum * vertNum];

	for (int i = 0; i < vertNum * vertNum; ++i) {
		mAdjacent[i] = 0;
	}

	//A sparse graph is generted by default edge number
	if (!edgeNum)
		edgeNum = vertNum;

	//force nonMultiEdge cnstraint to be false if there are more edges than a full graph
	if (nonMultiEdge && edgeNum > vertNum * vertNum)
		nonMultiEdge = false;

	for (int i = 0; i < edgeNum; ++i) {
		int u, v;
		do {
			u = rand() % vertNum;
			v = rand() % vertNum;
		}while (nonMultiEdge && mAdjacent[u + v * vertNum]);
		mAdjacent[u + v * vertNum] = 1;
		fout << u << ',' << v << endl;
	}
	fout.close();
}

