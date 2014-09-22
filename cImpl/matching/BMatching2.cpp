/*
 * BMatching2.cpp
 *
 *  Created on: 27/08/2014
 *      Author: jefcha
 */

#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <algorithm>
#include <cstdlib>
#include <cassert>



//
// DECLARATIONS
//

typedef std::vector<bool> BINARY_VECTOR;
typedef std::vector<int> INT_VECTOR;
typedef std::vector<float> FLOAT_VECTOR;
typedef std::vector<double> DOUBLE_VECTOR;

typedef std::list<int> INT_LIST;
typedef std::unordered_set<int> INT_SET;


void initAssignment(const FLOAT_VECTOR& vCostMat, int n, int m, INT_VECTOR& vAssignMat, INT_VECTOR& vRow, INT_VECTOR& vInvRow, FLOAT_VECTOR& vDualU, FLOAT_VECTOR& vDualV, INT_VECTOR& vUnassigneVerts);

float minCol(const FLOAT_VECTOR& vCostMat, int n, int m, int row);

float minRow(const FLOAT_VECTOR& vCostMat, int n, int m, int col);


int shortestPath(int currVert, const FLOAT_VECTOR& vCostMat, int n, int m,const INT_VECTOR& vRow, const INT_VECTOR& vInvRow, FLOAT_VECTOR& vDualU, FLOAT_VECTOR& vDualV, INT_VECTOR& vPred);


/* ****************************************************************** */

/**
 * Maximum matching algorithm for continuous costs.
 */
double matching2(int n, int m, const std::vector<float>& vA, std::vector<int>& m1, std::vector<int>& m2)
{
    using namespace std;

    double simVal = 0;


    INT_VECTOR vRow(m, 0);

    INT_VECTOR vInvRow(n, 0);

    // X assignment matrix
    INT_VECTOR vAssignMat(n*m, 0);

    // dual variables
    FLOAT_VECTOR vDualU(n, 0);
    FLOAT_VECTOR vDualV(m, 0);

    INT_VECTOR vUnassignedVerts;

    // initialise the variables (pass by reference)
    initAssignment(vA, n, m, vAssignMat, vRow, vInvRow, vDualU, vDualV, vUnassignedVerts);

    for (typename INT_VECTOR::const_iterator vit = vUnassignedVerts.begin(); vit != vUnassignedVerts.end(); ++vit) {
    	cout << *vit << endl;
    }


    for (typename INT_VECTOR::const_iterator vit = vUnassignedVerts.begin(); vit != vUnassignedVerts.end(); ++vit) {

    	int currVert = *vit;

    	// tree
    	INT_VECTOR vPred(n, 0);

    	int sinkVert = shortestPath(currVert, vA, n, m, vRow, vInvRow, vDualU, vDualV, vPred);
    	int j = sinkVert;

    	int i;
    	do {
    		i = vPred[j];
    		vRow[j] = i;
    		int h = vInvRow[i];
    		vInvRow[i] = j;
    		j = h;
    	} while(i != currVert);
    } // end of while loop

    // compute the sum of assigned costs

    return simVal;
} // end of function


/* ********************************************************* */

void initAssignment(const FLOAT_VECTOR& vCostMat, int n, int m, INT_VECTOR& vAssignMat, INT_VECTOR& vRow, INT_VECTOR& vInvRow, FLOAT_VECTOR& vDualU, FLOAT_VECTOR& vDualV, INT_VECTOR& vUnassignedVerts)
{

	for (int i = 0; i < n; ++i) {
		vDualU[i] = minCol(vCostMat, n, m, i);
		vUnassignedVerts.push_back(i);
	}

	// randomise the order
	std::random_shuffle(vUnassignedVerts.begin(), vUnassignedVerts.end());

	for (int j = 0; j < m; ++j) {
		vDualV[j] = minRow(vCostMat, n, m, j);
	}

	// assuming vAssignMat and vRow are all zeros
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			if (vRow[j] == 0 && (vCostMat[i + j*n] - vDualU[i] - vDualV[j] == 0)) {
				vAssignMat[i + j*n] = 1;
				vRow[j] = i;
				vInvRow[i] = j;
			}
		}
	}

} // end of initAssignment()



float minCol(const FLOAT_VECTOR& vCostMat, int n, int m, int row)
{
	assert(m > 0);

	float minVal = vCostMat[row + 0];

	for (int j = 1; j < m; ++j) {
		if (vCostMat[row + j * n] < minVal) {
			minVal = vCostMat[row + j * n];
		}
	}

	return minVal;
} // end of minCol




float minRow(const FLOAT_VECTOR& vCostMat, int n, int m, int col)
{
	assert(m > 0);

	float minVal = vCostMat[0 + col * n];

	for (int i = 1; i < n; ++i) {
		if (vCostMat[i + col * n] < minVal) {
			minVal = vCostMat[i + col * n];
		}
	}

	return minVal;
} // end of minRow


int shortestPath(int currVert, const FLOAT_VECTOR& vCostMat, int n, int m, const INT_VECTOR& vRow, const INT_VECTOR& vInvRow, FLOAT_VECTOR& vDualU, FLOAT_VECTOR& vDualV, INT_VECTOR& vPred)
{
	static const float large = 10^9;
	FLOAT_VECTOR vPi(m, large);

	INT_LIST vScannedRowVerts, vScannedColVerts;
	INT_SET vNotScannedRowVerts, vNotScannedColVerts;

	for (int i = 0; i < n; ++i) {
		vNotScannedRowVerts.insert(i);
	}

	for (int j = 0; j < m; ++j) {
		vNotScannedColVerts.insert(j);
	}

	int sinkVert = 0;
	float delta = 0;
	int iVert = currVert;

	while (sinkVert == 0) {
		vScannedRowVerts.push_back(iVert);
		vNotScannedRowVerts.erase(iVert);

		for (typename INT_SET::const_iterator jit = vNotScannedColVerts.begin(); jit != vNotScannedColVerts.end(); ++jit) {
			int jVert = *jit;
			if (delta + vCostMat[iVert + jVert * n] - vDualU[iVert] - vDualV[jVert] < vPi[jVert]) {
				vPred[jVert] = iVert;
				vPi[jVert] = delta + vCostMat[iVert + jVert * n] - vDualU[iVert] - vDualV[jVert];
			}
		}

		typename INT_SET::const_iterator hit = vNotScannedColVerts.begin();
		assert(hit != vNotScannedColVerts.end());
		int jVert = *hit;

		for ( ; hit != vNotScannedColVerts.end(); ++hit) {
			if (vPi[*hit] < vPi[jVert]) {
				jVert = *hit;
			}
		}

		vScannedColVerts.push_back(jVert);
		vNotScannedColVerts.erase(jVert);
		delta = vPi[jVert];

		if (vRow[jVert] == 0) {
			sinkVert = jVert;
		}
		else {
			iVert = vRow[jVert];
		}


	} // end of while loop

	// update the duals
	vDualU[currVert] += delta;
	for (typename INT_LIST::const_iterator uit = vScannedRowVerts.begin(); uit != vScannedRowVerts.end(); ++uit) {
		if (*uit != currVert) {
			vDualU[*uit] += delta - vPi[vInvRow[*uit]];
		}
	}

	for (typename INT_LIST::const_iterator vit = vScannedColVerts.begin(); vit != vScannedColVerts.end(); ++vit) {
		vDualV[*vit] +=  -delta + vPi[*vit];
	}

	return sinkVert;
} // end of augment()

