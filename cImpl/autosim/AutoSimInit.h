/*
 * RoleSim2Init.h
 *
 *  Created on: 13/08/2014
 *      Author: jefcha
 */

#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "../utils/pairUtils.h"

#ifndef ROLESIM2INIT_H_
#define ROLESIM2INIT_H_

/**
 * Abstract class for initialisation.
 */
class InitAlgor
{
protected:
	float m_dampingFactor;


public:
	InitAlgor(float dampingFactor);

	virtual ~InitAlgor();

	/**
	 * Initialisation for complete similarity matrix.
	 */
	virtual void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
			float* mSim, float ioBalFactor) = 0;

	/**
	 * Initialisation for complete similarity matrix but with initialisation for a bool matrix indicating which entries we need to compute.
	 */
	virtual void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
				float* mSim, float ioBalFactor, bool* mCompEntry) = 0;

	/**
	 * Initialisation for parital similarity matrix.
	 */
	virtual void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
			C_INT_PAIR_HASH_SET* hPairSim, C_INT_PAIR_HASH_MAP* vPrevSim, C_INT_PAIR_HASH_MAP* vCurrSim, float ioBalFactor) = 0;
}; // end of abstract class InitAlgor()


/* ******************************************************************** */

/**
 * Degree ratio similarity value intialisation.
 */
class DegRatioInitAlgor : public InitAlgor
{
public:
	DegRatioInitAlgor(float dampingFactor);

	virtual ~DegRatioInitAlgor();

	void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
			float* mSim, float ioBalFactor);

	void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
				float* mSim, float ioBalFactor, bool* mCompEntry);

	void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
			C_INT_PAIR_HASH_SET* hPairSim, C_INT_PAIR_HASH_MAP* vPrevSim, C_INT_PAIR_HASH_MAP* vCurrSim, float ioBalFactor);
};



/**
 * Binary similarity value intialisation.
 */
class BinaryInitAlgor : public InitAlgor
{
public:
	BinaryInitAlgor(float dampingFactor);

	virtual ~BinaryInitAlgor();

	void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
			float* mSim, float ioBalFactor);

	void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
				float* mSim, float ioBalFactor, bool* mCompEntry);

	void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
			C_INT_PAIR_HASH_SET* hPairSim, C_INT_PAIR_HASH_MAP* vPrevSim, C_INT_PAIR_HASH_MAP* vCurrSim, float ioBalFactor);

};


/**
 * Degree Binary similarity value intialisation.
 */
class DegBinaryInitAlgor : public InitAlgor
{
public:
	DegBinaryInitAlgor(float dampingFactor);

	virtual ~DegBinaryInitAlgor();

	void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
			float* mSim, float ioBalFactor);

	void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
					float* mSim, float ioBalFactor, bool* mCompEntry);

	void initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
			C_INT_PAIR_HASH_SET* hPairSim, C_INT_PAIR_HASH_MAP* vPrevSim, C_INT_PAIR_HASH_MAP* vCurrSim, float ioBalFactor);
};


#endif /* ROLESIM2INIT_H_ */

