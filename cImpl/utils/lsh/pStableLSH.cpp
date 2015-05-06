/*
 * pStableLSH.cpp
 *
 *  Created on: 06/05/2015
 *      Author: jefcha
 */

#include <iostream>
#include <chrono>
#include <random>
#include <vector>
#include <cassert>
#include <cmath>
#include "pStableLSH.h"


PStableLSH::PStableLSH(int pointDim, int projWidth, int projPerHashVal, int hashTableNum)
	: m_projWidth(projWidth), m_projPerHashVal(projPerHashVal), m_hashTableNum(hashTableNum), m_vHashFunc(NULL), m_vHashTable(NULL)
{
	// initialise the hash functions
	m_vHashFunc = new PStableHash[projPerHashVal];
	for (int k = 0; k < projPerHashVal; ++k) {
		m_vHashFunc[k] = PStableHash(pointDim, projWidth);
	}

	// initialise the hash table
	m_vHashTable = new HASH_TABLE[hashTableNum];
	for (int h = 0; h < hashTableNum; ++h) {
		m_vHashTable[h] = HASH_TABLE(projWidth);
	}

} // end of PStableLSH()



PStableLSH::~PStableLSH()
{
	delete[] m_vHashTable;
	delete[] m_vHashFunc;
} // end of ~PStableLSH()



void PStableLSH::insertPoint(const POINT& vPoint, int pointId)
{
	// TODO: assume one hash function and table for now
	int bucket = m_vHashFunc[0].hash(vPoint);
	m_vHashTable[0][bucket].push_back(pointId);
} // end of insertPoint()



const std::vector<int>& PStableLSH::getBucket(int hashTableIndex, int bucket) const
{
	// TODO: assume one hash table for now
	return m_vHashTable[0][bucket];
}







/* ************************************************************************** */


PStableHash::PStableHash(int pointDim, int projWidth)
	: m_projWidth(projWidth), m_vA(pointDim)
{
	// construct a trivial random generator engine from a time-based seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	// normal distribution with mean 0.0 and variance 1.0
	std::normal_distribution<double> gaussDist(0.0, 1.0);
	// compute the hash constants
	for (int i = 0; i < pointDim; ++i) {
		m_vA[i] = gaussDist(generator);
	}

	// initialise bias
	std::uniform_real_distribution<double> unifDist(0.0, projWidth);
	m_b = unifDist(generator);
} // end of PStableHash()


int PStableHash::hash(POINT vPoint) const
{
	assert(m_vA.size() == vPoint.size());

	// dot product
	double prod = 0.0;
	for (int i = 0; i < m_vA.size(); ++i) {
		prod += m_vA[i] * vPoint[i];
	}

	return std::floor((prod + m_b) / m_projWidth);
} // end of hash()
