/*
 * cRegeImplOpti.h
 *
 *  Created on: 05/10/2014
 *      Author: youhan
 */

#include <list>
#include <vector>

#include "../similarity/dampedSimilarity.h"

using namespace std;

#ifndef CREGEIMPLOPTI_H_
#define CREGEIMPLOPTI_H_

class RegeImplOpti : public IterSimilarity
{
protected:




public:

	RegeImplOpti(int maxIter);

	RegeImplOpti(int maxIter, float convEpsilon);

	virtual ~RegeImplOpti();

	virtual float* computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum);

	int index(int a, int b, int vertNum);

	void loopOverNeigh(int i, int j, float* mPrevSim, float* mCurrSim_i_j, int vertNum, vector< vector<int> > vvNeigh);

};



#endif /* CREGEIMPLOPTI_H_ */
