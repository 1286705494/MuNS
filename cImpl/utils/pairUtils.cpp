/*
 * pairUtils.cpp
 *
 *  Created on: 03/09/2014
 *      Author: jkct440
 */

#include <boost/functional/hash.hpp>
#include "pairUtils.h"


std::size_t PairHash::operator()(const std::pair<int, int>& vertPair) const
{
	std::size_t seed = 0;
	boost::hash_combine(seed, vertPair.first);
	boost::hash_combine(seed, vertPair.second);

	return seed;
} // end of operator()
