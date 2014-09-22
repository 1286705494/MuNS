/*
 * cAutoSimIceberg.h
 *
 *  Created on: 03/09/2014
 *      Author: jkct440
 */

#ifndef CAUTOSIMICEBERG_H_
#define CAUTOSIMICEBERG_H_

#include <list>
#include <string>
#include <unordered_set>
#include <stdexcept>

//#include "../similarity/dampedSimilarity.h"
#include "cAutoSim.h"
#include "AutoSimInit.h"
#include "../utils/pairUtils.h"



#ifdef _ZERO_SIM_COUNT
	extern
	std::list<int> G_vZeroSimCount;
#endif





//class PairEqual
//{
//public:
//
//};



class AutoSimIceberg : public AutoSim
{
protected:

	/** Iceberg filter parameter.  All similarities less than this are estimated. */
	float m_simThres;

	/** Iceberg approximation factor for non-computed similarities. */
	float m_approxFaction;


	/** Use early similarity stop. */
	bool m_bEarlySimStop2;
	float m_earlySimStopThres2;

	/** Set of vertices pairs that we need to compute similarities over. */
//	std::unordered_set<std::pair<int, int>, PairHash > m_hCompVertPairs;


public:

//	AutoSimIceberg(float dampingFactor, int maxIter, const std::string& sInitAlgor, float simThres, float approxFaction) throw(std::invalid_argument);
//	AutoSimIceberg(float dampingFactor, int maxIter, float convEpsilon, const std::string& sInitAlgor, float simThres, float approxFaction) throw(std::invalid_argument);

	AutoSimIceberg(float dampingFactor, int maxIter, const std::string& sInitAlgor, bool earlySimStop, float earlySimStopThres, bool useInputBalance, float ioBalance, float simThres, float approxFaction) throw(std::invalid_argument);
	AutoSimIceberg(float dampingFactor, int maxIter, float convEpsilon, const std::string& sInitAlgor, bool earlySimStop, float earlySimStopThres, bool useInputBalance, float ioBalance, float simThres, float approxFaction) throw(std::invalid_argument);

	virtual ~AutoSimIceberg();

	virtual float* computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum);




protected:

	/**
	 * Perform the iceberg filtering.
	 */
	void icebergFilter(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh, C_INT_PAIR_HASH_SET* hValidPair, float beta);


	/**
	 * Estimate similarity value.
	 */
	float estimateSim(int vert1, int vert2, const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh, float beta) const;

	/**
	 * Filter out pairs whose similarity become less than m_earlySimStopThres2.
	 */
	void filterEarlySimPairs(C_INT_PAIR_HASH_SET* phValidPair, const C_INT_PAIR_HASH_MAP* const pmPrevSim, const C_INT_PAIR_HASH_MAP* const pmCurrSim);


	/**
	 * LI distance between two matrices represented as hash tables.
	 */
	float matDiff(const C_INT_PAIR_HASH_SET* const phValidPair, const C_INT_PAIR_HASH_MAP* const pmPrevSim, const C_INT_PAIR_HASH_MAP* const pmCurrSim) const;

}; // end of class AutoSimIceberg







#endif /* CAUTOSIMICEBERG_H_ */
