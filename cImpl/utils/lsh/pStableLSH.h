/*
 * pStableLSH.h
 *
 * Locality sensitive hashing for euclidean distance, using p-stable distribution as the hashing family.
 *
 *  Created on: 06/05/2015
 *      Author: jefcha
 */

#ifndef PSTABLELSH_H_
#define PSTABLELSH_H_

#include <vector>
#include "LSH.h"




// advanced declaration
class PStableHash;



/**
 * LSH based on the P-Stable distribution for L2 norm.
 */
class PStableLSH : public LSH
{
protected:

typedef std::vector<std::vector<int> > HASH_TABLE;

protected:

	/* Parameters. */

	/** Width of projection (r). */
	int m_projWidth;

	/** Projections per hash value (k). */
	int m_projPerHashVal;

	/** Number of hash tables (l). */
	int m_hashTableNum;


	/** Pointer of Hash functions. */
	PStableHash** m_vHashFunc;

	/** Hash tables. */
	HASH_TABLE* m_vHashTable;




public:

	PStableLSH(int pointDim, int projWidth, int projPerHashVal = 1, int hashTableNum = 1);

	~PStableLSH();

	/**
	 * Insert a point into the hash table.
	 */
	virtual void insertPoint(const POINT& vPoint, int pointId);


	/**
	 * Get the contents of bucket 'bucket' from hash table 'hashTableIndex'.
	 */
	const std::vector<int>& getBucket(int hashTableIndex, int bucket) const;



protected:



}; // end of class PStableLSH


/* ************************************************************************** */


/**
 * Hash function based on P-Stable distribution.
 */
class PStableHash
{
protected:
	/** Projection width (r). */
	int m_projWidth;

	/** Internal 'a' vector. */
	std::vector<double> m_vA;

	/** Offset value. */
	double m_b;

public:

	PStableHash(int pointDim, int projWidth);

	int hash(POINT vPoint) const;

}; // end of class PStableHash


#endif /* PSTABLELSH_H_ */
