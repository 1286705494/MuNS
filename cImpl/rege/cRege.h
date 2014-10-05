/*
 * cRege.h
 *
 *  Created on: 05/10/2014
 *      Author: youhan
 */

#include <list>
#include <vector>

#include "../similarity/dampedSimilarity.h"

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


};



#endif /* CREGE_H_ */
