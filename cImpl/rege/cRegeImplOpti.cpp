/*
 * cRegeImplOpti.cpp
 *
 *  Created on: 05/10/2014
 *      Author: youhan
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <cassert>
#include <cmath>

#include "cRegeImplOpti.h"
#include "../utils/matUtils.h"

#define ITER typename vector<int>::const_iterator

using namespace std;

//inline int RegeImplOpti::index(int r, int c, int vertNum) { return r < c ? r + c * vertNum : c + r * vertNum; }

inline int RegeImplOpti::index(int a, int b, int vertNum)
{
	int sum = a + b;
	int diff = abs(a - b);
	int c2 = (sum + diff);
	int r2 = (sum - diff);
	return ((c2 * (c2 + 2)) / 4 + r2) / 2;
}

void RegeImplOpti::loopOverNeigh(int i, int j, float* mPrevSim, float* mCurrSim_i_j, int vertNum, vector< vector<int> > vvNeigh)
{

	if (vvNeigh[i].size() == 0 || vvNeigh[j].size() == 0) return;
	for (ITER kit = vvNeigh[i].begin(); kit != vvNeigh[i].end(); ++kit) {
		ITER lit = vvNeigh[j].begin();
		int lmax = index(*kit, *lit, vertNum);
		for (; lit != vvNeigh[j].end(); ++lit) {
			if (mPrevSim[lmax] < mPrevSim[index(*kit, *lit, vertNum)])
				lmax = index(*kit, *lit, vertNum);
		} // end of inner for
		*mCurrSim_i_j += mPrevSim[lmax];
	} // end of outer for

}

RegeImplOpti::RegeImplOpti(int maxIter)
	: IterSimilarity(maxIter)
{
} // end of Rege()


RegeImplOpti::RegeImplOpti(int maxIter, float convEpsilon)
	: IterSimilarity(maxIter, convEpsilon)
{
} // end of Rege()

RegeImplOpti::~RegeImplOpti()
{
}


float* RegeImplOpti::computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum)
{

	int mSize = (vertNum + 1) * vertNum / 2;

	// similarity matrix (column-major, and we only store the upper triangle, minus the diagonal (self-similarities, which is always 1)
	float* mPrevSim = new float[mSize];
	float* mCurrSim = new float[mSize];

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

	// initialise the values of simRank
	// all to 1
	for (int i = 0; i < mSize; ++i) {
			// no need to use index() as r always < c)
			mPrevSim[i] = 1;
	}

	// diagonals to 1 (since this isn't updated afterwards)
	for (int i = 0; i < vertNum; ++i) {
		mCurrSim[index(i, i, vertNum)]  = 1;
	}

	// perform loop iterations
	float* mTempSim;

	for (int t = 1; t <= m_maxIter; ++t) {
		cout << "iteration " << t << endl;

		// loop through all upper non-diagonal pairs
		for (int i = 0; i < vertNum; ++i) {
			for (int j = i+1; j < vertNum; ++j) {
				int index_i_j = index(i, j, vertNum);

				mCurrSim[index_i_j] = 0;

				// loop over the neighbours
				loopOverNeigh(i, j, mPrevSim, mCurrSim + index_i_j, vertNum, vvInNeigh);
				loopOverNeigh(i, j, mPrevSim, mCurrSim + index_i_j, vertNum, vvOutNeigh);
				loopOverNeigh(j, i, mPrevSim, mCurrSim + index_i_j, vertNum, vvInNeigh);
				loopOverNeigh(j, i, mPrevSim, mCurrSim + index_i_j, vertNum, vvOutNeigh);


				mCurrSim[index_i_j] /= vvInNeigh[i].size() + vvOutNeigh[i].size() + vvInNeigh[j].size() + vvOutNeigh[j].size();

			} // end of inner for
		} // end of outer for

		// loop through diagonal pairs (set to one) (no need to, since it doesn't change)


#ifdef _COLLECT_SIM_DELTA_
		// do comparison
		if (t > 1) {
			assert(mPrevSim != NULL);
			assert(mCurrSim != NULL);

			m_vSimDelta.push_back(l1MatDiff(mPrevSim, mCurrSim, vertNum, vertNum));
		}
#endif

		m_iterRan = t;

		// check if we have convergence epsilon
		if (m_bUseConvEpsilon) {
//			cout << "l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum) = " << l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum) << endl;
			if (l1MatDiff(mPrevSim, mCurrSim, vertNum, vertNum) / (vertNum * vertNum) < m_convEpsilon) {
				break;
			}
		}

		// swap the sim matrices so we do not need to allocate more matrices then need to be
		mTempSim = mCurrSim;
		mCurrSim = mPrevSim;
		mPrevSim = mTempSim;

	} // end of loop through iterations
	// destroy dynamically allocated memory

	delete[] mCurrSim;

	float * mReturn = new float[vertNum * vertNum];

	for (int i = 0; i < vertNum; ++i) {
		for (int j = 0; j < vertNum; ++j) {
			mReturn[i * vertNum + j] = mPrevSim[index(i, j, vertNum)];
		}
	}

	delete[] mPrevSim;

	return mReturn;
}
