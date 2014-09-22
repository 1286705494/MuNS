/*
 * pairUtils.h
 *
 *  Created on: 03/09/2014
 *      Author: jkct440
 */

#include <unordered_set>
#include <unordered_map>

#ifndef PAIRUTILS_H_
#define PAIRUTILS_H_

class PairHash
{
public:
	size_t operator() (const std::pair<int, int>& vertPair) const;
};


// Some common pairs datastructures
typedef std::unordered_map<std::pair<int, int>, float, PairHash> C_INT_PAIR_HASH_MAP;
typedef std::unordered_set<std::pair<int, int>, PairHash> C_INT_PAIR_HASH_SET;


#endif /* PAIRUTILS_H_ */

