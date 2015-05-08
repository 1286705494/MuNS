/*
 * cAutoSimHash.h
 *
 *  Created on: 19/04/2015
 *      Author: jkct440
 */

#ifndef CAUTOSIMHASH_H_
#define CAUTOSIMHASH_H_


#include <list>
#include <string>
#include <unordered_set>
#include <stdexcept>

//#include "../similarity/dampedSimilarity.h"
#include "cAutoSim.h"
#include "AutoSimInit.h"
#include "../utils/pairUtils.h"
#include "../utils/lsh/LSH.h"
#include "../utils/lsh/pStableLSH.h"



#ifdef _ZERO_SIM_COUNT
	extern
	std::list<int> G_vZeroSimCount;
#endif





class AutoSimHash : public AutoSim
{
protected:

	/** Iceberg approximation factor for non-computed similarities. */
	float m_approxFaction;


	/** Use early similarity stop. */
	bool m_bEarlySimStop2;
	float m_earlySimStopThres2;

	/** Hashing parameters. */
	// Number of bins/length to represent the neighbourhood degree vectors
	int m_binNum;
	// Number of hash buckets
	int m_hashBucketNum;

	/** Hashing function. */
	LSH* m_hashFunc;




public:

	AutoSimHash(float dampingFactor, int maxIter, const std::string& sInitAlgor, bool earlySimStop, float earlySimStopThres, bool useInputBalance, float ioBalance, float approxFaction, const std::string& sHashFunction, int binNum, int hashBucketNum) throw(std::invalid_argument);
	AutoSimHash(float dampingFactor, int maxIter, float convEpsilon, const std::string& sInitAlgor, bool earlySimStop, float earlySimStopThres, bool useInputBalance, float ioBalance, float approxFaction, const std::string& sHashFunction, int binNum, int hashBucketNum) throw(std::invalid_argument);

	virtual ~AutoSimHash();

	virtual float* computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum);


protected:

	/**
	 * Perform the iceberg filtering.
	 */
	void hashFilter(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh, C_INT_PAIR_HASH_SET* hValidPair, int maxInDeg, int maxOutDeg);

	/**
	 * Hash.
	 */
	void hash(std::vector<std::vector<int> >& vvVertBuckets, const std::vector<int>& vInHist, const std::vector<int>& vOutHist);


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


	/**
	 * Construct degree vector.
	 */
	void constructNeighDegVec(std::vector<std::vector<int> >& vvInHist, std::vector<std::vector<int> >& vvOutHist, const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh, int maxInDeg, int maxOutDeg) const;


	/**
	 * LSH Hash the neighbourhood degree vectors into bins.
	 */
	void hashVectors(const std::vector<std::vector<int> >& vvInHist, const std::vector<std::vector<int> >& vvOutHist) const;

}; // end of class AutoSimHash



#endif /* CAUTOSIMHASH_H_ */
