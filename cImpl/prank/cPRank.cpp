/*
 * cPRank.cpp
 *
 *  Created on: 17/08/2014
 *      Author: jkct440
 */

#include <iostream>
#include <vector>
#include <list>
#include <cassert>

#include "cPRank.h"
#include "../utils/matUtils.h"




PRank::PRank(float dampingFactor, int maxIter, float ioBalance)
	: DampedSimilarity(dampingFactor, maxIter), m_ioBalance(0.5)
{
} // end of PRank()


PRank::PRank(float dampingFactor, int maxIter, float convEpsilon, float ioBalance)
	: DampedSimilarity(dampingFactor, maxIter, convEpsilon), m_ioBalance(0.5)
{
} // end of PRank()


PRank::~PRank()
{
}


float* PRank::computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum)
{
	using namespace std;

    // similarity matrix (column-major)
    float* mPrevSim = new float[vertNum*vertNum];
    float* mCurrSim = new float[vertNum*vertNum];
    float **pmPrevSim = &mPrevSim;
    float **pmCurrSim = &mCurrSim;

    // construct neighbour list
    vector< vector<int> > vvInNeigh(vertNum);
    vector< vector<int> > vvOutNeigh(vertNum);

    // set the neighbourhoods and degrees
    std::list<int>::const_iterator sit = vSrc.begin(), tit = vTar.begin();
    for ( ; sit != vSrc.end(); ++sit, ++tit) {
    	assert(*sit < vertNum && *tit < vertNum);
        vvInNeigh[*tit].push_back(*sit);
        vvOutNeigh[*sit].push_back(*tit);
    } // end of for


    // initialise the values of pRank
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

	    		// loop over the in neighbours
	    		typename vector<int>::const_iterator niti = vvInNeigh[i].begin();
	    		typename vector<int>::const_iterator nitj;

	    		float simTotal = 0;
	    		for ( ; niti != vvInNeigh[i].end(); ++niti) {
	    			nitj = vvInNeigh[j].begin();
	    			for ( ; nitj != vvInNeigh[j].end(); ++nitj) {
	    				simTotal += mTempPrevSim[*niti + *nitj * vertNum];
	    			}
	    		}

	    		if (vvInNeigh[i].size() > 0 && vvInNeigh[j].size() > 0) {
	    			mTempCurrSim[i + j * vertNum] = m_ioBalance * m_dampingFactor * simTotal / (vvInNeigh[i].size() * vvInNeigh[j].size());
	    		}
	    		else {
	    			mTempCurrSim[i + j * vertNum] = 0;
	    		}

	    		// loop over the out neighbours
	    		niti = vvOutNeigh[i].begin();

	    		simTotal = 0;
	    		for ( ; niti != vvOutNeigh[i].end(); ++niti) {
	    			nitj = vvOutNeigh[j].begin();
	    			for ( ; nitj != vvOutNeigh[j].end(); ++nitj) {
	    			    simTotal += mTempPrevSim[*niti + *nitj * vertNum];
	    			}
	    		}

	    		if (vvOutNeigh[i].size() > 0 && vvOutNeigh[j].size() > 0) {
	    			mTempCurrSim[i + j * vertNum] += (1 - m_ioBalance) * m_dampingFactor * simTotal / (vvOutNeigh[i].size() * vvOutNeigh[j].size());
	    		}
	    		else {
	    			mTempCurrSim[i + j * vertNum] += 0;
	    		}
	    		mTempCurrSim[j + i * vertNum] = mTempCurrSim[i + j * vertNum];
            } // end of inner for
        } // end of outer for


	    // loop through diagonal pairs (always 1)
	    for (int i = 0; i < vertNum; ++i) {

	    	// loop over in neighbours
//	    	typename vector<int>::const_iterator niti1 = vvInNeigh[i].begin();
//	    	typename vector<int>::const_iterator niti2 = vvInNeigh[i].begin();
//
//	    	float simTotal = 0;
//	    	for ( ; niti1 != vvInNeigh[i].end(); ++niti1) {
//	    		for ( ; niti2 != vvInNeigh[i].end(); ++niti2) {
//	    			simTotal += mTempPrevSim[*niti1 + *niti2 * vertNum];
//	    		}
//    		}
//
//    		mTempCurrSim[i + i * vertNum] = m_ioBalance * m_dampingFactor * simTotal / (vvInNeigh[i].size() * vvInNeigh[i].size());
//
//    		// loop over out neighbours
//	    	niti1 = vvOutNeigh[i].begin();
//	    	niti2 = vvOutNeigh[i].begin();
//
//	    	simTotal = 0;
//	    	for ( ; niti1 != vvOutNeigh[i].end(); ++niti1) {
//	    		for ( ; niti2 != vvOutNeigh[i].end(); ++niti2) {
//	    			simTotal += mTempPrevSim[*niti1 + *niti2 * vertNum];
//	    		}
//    		}

    		mTempCurrSim[i + i * vertNum] = 1;
        } // end of outer for





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

