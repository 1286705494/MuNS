/*
 * cSimRank.cpp
 *
 *  Created on: 16/08/2014
 *      Author: jefcha
 */

#include <iostream>
#include <vector>
#include <list>
#include <cassert>

#include "cSimRank.h"
#include "../utils/matUtils.h"


SimRank::SimRank(float dampingFactor,int maxIter)
	: DampedSimilarity(dampingFactor, maxIter)
{
} // end of SimRank()


SimRank::SimRank(float dampingFactor,int maxIter, float convEpsilon)
	: DampedSimilarity(dampingFactor, maxIter, convEpsilon)
{
} // end of SimRank()


SimRank::~SimRank()
{
}


float* SimRank::computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum)
{
	using namespace std;

    // similarity matrix (column-major)
	float* mPrevSim = new float[vertNum*vertNum];
	float* mCurrSim = new float[vertNum*vertNum];
	float **pmPrevSim = &mPrevSim;
	float **pmCurrSim = &mCurrSim;

    // construct neighbour list
    vector< vector<int> > vvInNeigh(vertNum);

    // set the neighbourhoods and degrees
    std::list<int>::const_iterator sit = vSrc.begin(), tit = vTar.begin();
    for ( ; sit != vSrc.end(); ++sit, ++tit) {
    	assert(*sit < vertNum && *tit < vertNum);
        vvInNeigh[*tit].push_back(*sit);
    } // end of for


    // initialise the values of simRank
    // non-diagonals to 0
    for (int c = 0; c < vertNum; ++c) {
    	for (int r = 0; r < vertNum; ++r) {
    		mPrevSim[r + c*vertNum] = 0;
    	}
    }
    // diagonals to 1
    for (int i = 0; i < vertNum; ++i) {
    	mPrevSim[i + i * vertNum]  = 1;
    }

    // perform loop iterations
    for (int t = 1; t <= m_maxIter; ++t) {
    	cout << "iteration " << t << endl;

    	float* mTempPrevSim = *pmPrevSim;
    	float *mTempCurrSim = *pmCurrSim;

        // loop through all non-diagonal pairs
	    for (int i = 0; i < vertNum; ++i) {
	    	for (int j = i+1; j < vertNum; ++j) {

	    		// loop over the neighbours
	    		typename vector<int>::const_iterator niti = vvInNeigh[i].begin();
	    		typename vector<int>::const_iterator nitj = vvInNeigh[j].begin();

	    		float simTotal = 0;
	    		for ( ; niti != vvInNeigh[i].end(); ++niti) {
	    			for ( ; nitj != vvInNeigh[j].end(); ++nitj) {
	    				simTotal += mTempPrevSim[*niti + *nitj * vertNum];
	    			}
	    		}

	    		if (vvInNeigh[i].size() > 0 && vvInNeigh[j].size() > 0) {
	    			mTempCurrSim[i + j * vertNum] = m_dampingFactor * simTotal / (vvInNeigh[i].size() * vvInNeigh[j].size());
	    		}
	    		else {
	    			mTempCurrSim[i + j * vertNum] = 0;
	    		}
	    		mTempCurrSim[j + i * vertNum] = mTempCurrSim[i + j * vertNum];
            } // end of inner for
        } // end of outer for

	    // loop through diagonal pairs (set to one) (no need to, since it doesn't change
	    for (int i = 0; i < vertNum; ++i) {
	    	mTempCurrSim[i + i * vertNum] = 1;
	    }


#ifdef _COLLECT_SIM_DELTA_
    	// do comparison
        if (t > 1) {
        	assert(*pmPrevSim != NULL);
        	assert(*pmCurrSim != NULL);

        	m_vSimDelta.push_back(l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum));
        }
#endif

        m_iterRan = t;

        // check if we have convergence epsilon
        if (m_bUseConvEpsilon) {
//        	cout << "l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum) = " << l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum) << endl;
        	if (l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum) / (vertNum * vertNum) < m_convEpsilon) {
        		break;
        	}
        }



        // swap the sim matrices so we do not need to allocate more matrices then need to be
        float **pmTempSim = pmPrevSim;
        pmPrevSim = pmCurrSim;
        pmCurrSim = pmTempSim;


    } // end of loop through iterations

    // destroy dynamically allocated memory
    delete[] *pmPrevSim;
    pmPrevSim = NULL;

    return *pmCurrSim;
} // end of computeSim()


