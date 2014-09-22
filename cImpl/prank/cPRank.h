/*
 * cPRank.h
 *
 *  Created on: 17/08/2014
 *      Author: jkct440
 */

#include "../similarity/dampedSimilarity.h"

#ifndef CPRANK_H_
#define CPRANK_H_



class PRank : public DampedSimilarity
{
protected:


	float m_ioBalance;




public:
	PRank(float dampingFactor, int maxIter, float ioBalance);

	PRank(float dampingFactor, int maxIter, float convEpsilon, float ioBalance);

	virtual ~PRank();

	virtual float* computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum);


};



#endif /* CPRANK_H_ */
