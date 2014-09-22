/*
 * cSimRank.h
 *
 *  Created on: 16/08/2014
 *      Author: jefcha
 */

#include <list>
#include <vector>

#include "../similarity/dampedSimilarity.h"

#ifndef CSIMRANK_H_
#define CSIMRANK_H_

class SimRank : public DampedSimilarity
{
protected:




public:

	SimRank(float dampingFactor, int maxIter);

	SimRank(float dampingFactor, int maxIter, float convEpsilon);

	virtual ~SimRank();

	virtual float* computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum);


};




#endif /* CSIMRANK_H_ */
