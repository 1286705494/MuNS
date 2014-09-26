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

inline int index(int r, int c, int vertNum) { return r < c ? r + c * vertNum : c + r * vertNum; }



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

    // similarity matrix (column-major, and we only store the upper triangle, minus the diagonal (self-similarities, which is always 1)
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
    for (int r = 0; r < vertNum; ++r) {
    	for (int c = r+1; c < vertNum; ++c) {
    		// no need to use index() as r always < c)
    		mPrevSim[r + c*vertNum] = 0;
    	}
    }

    // diagonals to 1 (since this isn't updated afterwards)
    for (int i = 0; i < vertNum; ++i) {
    	mPrevSim[i + i * vertNum]  = 1;
    	mCurrSim[i + i * vertNum]  = 1;
    }

    // perform loop iterations
    for (int t = 1; t <= m_maxIter; ++t) {
    	cout << "iteration " << t << endl;

    	float* mTempPrevSim = *pmPrevSim;
    	float *mTempCurrSim = *pmCurrSim;

        // loop through all upper non-diagonal pairs
	    for (int i = 0; i < vertNum; ++i) {
	    	for (int j = i+1; j < vertNum; ++j) {

	    		// loop over the neighbours
	    		typename vector<int>::const_iterator niti = vvInNeigh[i].begin();

	    		float simTotal = 0;
	    		for ( ; niti != vvInNeigh[i].end(); ++niti) {
	    			typename vector<int>::const_iterator nitj = vvInNeigh[j].begin();
	    			for ( ; nitj != vvInNeigh[j].end(); ++nitj) {
	    				simTotal += mTempPrevSim[*niti + *nitj * vertNum];
//	    				simTotal += mTempPrevSim[index(*niti, *nitj, vertNum)];
	    			}
	    		}

	    		if (vvInNeigh[i].size() > 0 && vvInNeigh[j].size() > 0) {
	    			// no need to use index, as i < j
	    			mTempCurrSim[i + j * vertNum] = m_dampingFactor * simTotal / (vvInNeigh[i].size() * vvInNeigh[j].size());
	    		}
	    		else {
	    			// no need to use index, as i < j
	    			mTempCurrSim[i + j * vertNum] = 0;
	    		}
	    		mTempCurrSim[j + i * vertNum] = mTempCurrSim[i + j * vertNum];
            } // end of inner for
        } // end of outer for

	    // loop through diagonal pairs (set to one) (no need to, since it doesn't change)


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


