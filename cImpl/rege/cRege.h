/*
 * cRege.h
 *
 *  Created on: 05/10/2014
 *      Author: youhan
 */

#include <list>
#include <vector>

#include "../similarity/dampedSimilarity.h"

using namespace std;

#ifndef CREGE_H_
#define CREGE_H_

class Rege : public IterSimilarity
{
protected:




public:

	Rege(int maxIter);

	Rege(int maxIter, float convEpsilon);

	virtual ~Rege();

	virtual float* computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum);

	void loopOverNeigh(int i, int j, float* mPrevSim, float* mCurrSim, int vertNum, vector< vector<int> > vvNeigh);

};



#endif /* CREGE_H_ */
