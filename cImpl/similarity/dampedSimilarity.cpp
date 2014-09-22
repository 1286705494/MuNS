/*
 * dampedSimilarity.cpp
 *
 *  Created on: 18/08/2014
 *      Author: jefcha
 */


#include <vector>
#include <cassert>
#include "dampedSimilarity.h"









/* **************************************************** */

IterSimilarity::IterSimilarity(int maxIter)
	: m_maxIter(maxIter), m_convEpsilon(0.1), m_bUseConvEpsilon(false), m_iterRan(0)
{
}

IterSimilarity::IterSimilarity(int maxIter, float convEpsilon)
	: m_maxIter(maxIter), m_convEpsilon(convEpsilon), m_bUseConvEpsilon(true), m_iterRan(0)
{
}






IterSimilarity::~IterSimilarity()
{
}


void IterSimilarity::setMaxIter(int maxIter) {
	m_maxIter = maxIter;
}


int IterSimilarity::getIterRan() const
{
	return m_iterRan;
}


#ifdef _COLLECT_SIM_DELTA_

const std::vector<float>& IterSimilarity::getSimDelta() const
{
	return m_vSimDelta;
}

#endif


/* **************************************************** */

DampedSimilarity::DampedSimilarity(float dampedFactor, int maxIter)
	: IterSimilarity(maxIter), m_dampingFactor(dampedFactor)
{
	assert(m_dampingFactor >= 0 && m_dampingFactor <= 1);
}


DampedSimilarity::DampedSimilarity(float dampedFactor, int maxIter, float convEpsilon)
	: IterSimilarity(maxIter, convEpsilon), m_dampingFactor(dampedFactor)
{
	assert(m_dampingFactor >= 0 && m_dampingFactor <= 1);
}


DampedSimilarity::~DampedSimilarity()
{
}

void DampedSimilarity::setDampingFactor(float dampingFactor) {
	m_dampingFactor = dampingFactor;
}




