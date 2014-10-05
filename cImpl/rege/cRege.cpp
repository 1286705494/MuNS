/*
 * cRege.cpp
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

#include "cRege.h"
#include "../utils/matUtils.h"

#define ITER typename vector<int>::const_iterator

using namespace std;

//inline int index(int r, int c, int vertNum) { return r < c ? r + c * vertNum : c + r * vertNum; }

inline int index(int a, int b, int vertNum)
{
	int sum = a + b;
	int diff = abs(a - b);
	int c2 = (sum + diff);
	int r2 = (sum - diff);
	return ((c2 * (c2 + 2)) / 2 + r2) / 2;
}

void loopOverNeigh(int i, int j, float* mPrevSim, float* mCurrSim, int vertNum, vector< vector<int> > vvNeigh)
{
	if (vvNeigh[i].size() == 0 || vvNeigh[j].size() == 0) return;
	for (ITER kit = vvNeigh[i].begin(); kit != vvNeigh[i].end(); ++kit) {
		ITER lmax = vvNeigh[j].begin();
		for (ITER lit = vvNeigh[j].begin(); lit != vvNeigh[j].end(); ++lit) {
			if (mPrevSim[*kit + *lmax * vertNum] < mPrevSim[*kit + *lit * vertNum])
				lmax = lit;
		} // end of inner for
		mCurrSim[i + j * vertNum] += mPrevSim[*kit + *lmax * vertNum];
//		cout << i << "->" << *kit << " == " << j << "->" << *lmax << " in " << mPrevSim[*kit + *lmax * vertNum] << endl;
//		cout << "support[" << i << "][" << j <<"] == " << mCurrSim << endl;
//		cin.ignore();
	} // end of outer for

}

Rege::Rege(int maxIter)
	: IterSimilarity(maxIter)
{
} // end of Rege()


Rege::Rege(int maxIter, float convEpsilon)
	: IterSimilarity(maxIter, convEpsilon)
{
} // end of Rege()

Rege::~Rege()
{
}


float* Rege::computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum)
{


	// similarity matrix (column-major, and we only store the upper triangle, minus the diagonal (self-similarities, which is always 1)
	float* mPrevSim = new float[vertNum*vertNum];
	float* mCurrSim = new float[vertNum*vertNum];

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
	for (int r = 0; r < vertNum; ++r) {
		for (int c = 0; c < vertNum; ++c) {
			// no need to use index() as r always < c)
			mPrevSim[r + c*vertNum] = 1;
		}
	}

	// diagonals to 1 (since this isn't updated afterwards)
	for (int i = 0; i < vertNum; ++i) {
		mCurrSim[i + i * vertNum]  = 1;
	}

	// perform loop iterations
	float* mTempSim;

	for (int t = 1; t <= m_maxIter; ++t) {
		cout << "iteration " << t << endl;

		// loop through all upper non-diagonal pairs
		for (int i = 0; i < vertNum; ++i) {
			for (int j = i+1; j < vertNum; ++j) {

				mCurrSim[i + j * vertNum] = 0;
				mCurrSim[j + i * vertNum] = 0;

				// loop over the neighbours
				loopOverNeigh(i, j, mPrevSim, mCurrSim, vertNum, vvInNeigh);
				loopOverNeigh(i, j, mPrevSim, mCurrSim, vertNum, vvOutNeigh);
				loopOverNeigh(j, i, mPrevSim, mCurrSim, vertNum, vvInNeigh);
				loopOverNeigh(j, i, mPrevSim, mCurrSim, vertNum, vvOutNeigh);

//				cout << "m[" << i << "][" <<  j << "] = " << mCurrSim[i + j * vertNum] <<" / " << vvInNeigh[i].size() + vvOutNeigh[i].size() + vvInNeigh[j].size() + vvOutNeigh[j].size() << endl;
//				cin.ignore();

				mCurrSim[i + j * vertNum] += mCurrSim[j + i * vertNum];
				mCurrSim[i + j * vertNum] /= vvInNeigh[i].size() + vvOutNeigh[i].size() + vvInNeigh[j].size() + vvOutNeigh[j].size();
				mCurrSim[j + i * vertNum] = mCurrSim[i + j * vertNum];

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

//		for (int i = 0; i < vertNum ; ++i){
//			for (int j = 0; j < vertNum; ++j){
//				cout << fixed << setw(5) << setprecision(3) << mCurrSim[i + j * vertNum] << ' ';
//			}
//			cout <<endl;
//		}


		// swap the sim matrices so we do not need to allocate more matrices then need to be
		mTempSim = mCurrSim;
		mCurrSim = mPrevSim;
		mPrevSim = mTempSim;

	} // end of loop through iterations
	// destroy dynamically allocated memory
	delete[] mCurrSim;

	return mPrevSim;
}
