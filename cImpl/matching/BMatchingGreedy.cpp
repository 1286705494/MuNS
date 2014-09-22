/*
 * BMatchingGreedy.cpp
 *
 *  Created on: 07/08/2014
 *      Author: jefcha
 */

#include <vector>
#include <deque>
#include <list>
#include <cassert>
#include <algorithm>
#include "BMatchingGreedy.h"

float greedyMatching(int n, int m, const std::vector<float>& vA, int* m1, int* m2)
{
	using namespace std;

	typedef list<int> C_VertexList;
	typedef deque<int> C_VertexQueue;

	C_VertexQueue vRowUnmatched;
	C_VertexList vColUnmatched;
	for (int i = 0; i < n; ++i) {
		vRowUnmatched.push_back(i);
		vColUnmatched.push_back(i);
	}

	// random shuffle the rows
	random_shuffle(vRowUnmatched.begin(), vRowUnmatched.end());

	float simVal = 0;

	// loop through the rows (which has been randomised)
	int r = 0;
	while (vRowUnmatched.size() > 0) {
		int currRow = vRowUnmatched.front();
		vRowUnmatched.pop_front();

		typename C_VertexList::iterator maxColIt = vColUnmatched.end();
		float maxVal = -1;
		// find the maximum edge to one of nodes in vColUnmatched
		for (typename C_VertexList::iterator cit = vColUnmatched.begin(); cit != vColUnmatched.end(); ++cit) {
			float currVal = vA[currRow + *cit * n];
			if (currVal > maxVal) {
				maxColIt = cit;
				maxVal = currVal;
			}
		}

		simVal += maxVal;
		// make assignment
		m1[r] = currRow;
		m2[r] = *maxColIt;
		++r;

		// delete the matched column element
		vColUnmatched.erase(maxColIt);
	}

	return simVal;
} // end of greedyMatching()


