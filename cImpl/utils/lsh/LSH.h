/*
 * LSH.h
 *
 *  Created on: 08/05/2015
 *      Author: jefcha
 */

#ifndef LSH_H_
#define LSH_H_


/** Point type. */
typedef std::vector<int> POINT;
typedef std::vector<std::vector<int> > HASH_TABLE;

/**
 * Virtual hashing class.
 */
class LSH
{
public:

	virtual void insertPoint(const POINT& vPoint, int pointId) = 0;

	/**
	 * @returns: combined hash table.
	 */
	virtual const HASH_TABLE& getBuckets() const = 0;





}; // end of LSH



#endif /* LSH_H_ */
