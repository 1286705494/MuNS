/*
 * RoleSim2Init.cpp
 *
 *  Created on: 13/08/2014
 *      Author: jefcha
 */

#include <algorithm>
#include <utility>
#include <vector>
#include <unordered_map>
#include <cstdlib>

#include "AutoSimInit.h"
#include "../utils/pairUtils.h"

InitAlgor::InitAlgor(float dampingFactor)
	: m_dampingFactor(dampingFactor)
{
}

InitAlgor::~InitAlgor()
{
}

/* ***************************************************************** */

DegRatioInitAlgor::DegRatioInitAlgor(float dampingFactor)
	: InitAlgor(dampingFactor)
{
}

DegRatioInitAlgor::~DegRatioInitAlgor()
{}


void DegRatioInitAlgor::initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
		float* mSim,  float ioBalFactor)
{
	using namespace std;

	int vertNum = vvInNeigh.size();

    for (int i = 0; i < vertNum; ++i) {
        for (int j = 0; j < vertNum; ++j) {
        	int inSizei = vvInNeigh[i].size();
        	int inSizej = vvInNeigh[j].size();
        	int outSizei = vvOutNeigh[i].size();
        	int outSizej = vvOutNeigh[j].size();
            int maxInDeg = max(inSizei, inSizej);
            int maxOutDeg = max(outSizei, outSizej);
            float minInDeg = min(inSizei, inSizej);
            float minOutDeg = min(outSizei, outSizej);

            if (maxInDeg == 0) {
                if (maxOutDeg == 0) {
                	mSim[i + j*vertNum] = 0;
                }
                else {
                	float OUT = minOutDeg / maxOutDeg;
                    mSim[i + j*vertNum] = m_dampingFactor * OUT + 1 - m_dampingFactor;
                }
            }
            else {
            	float IN = minInDeg / maxInDeg;
                if (maxOutDeg == 0) {
                	mSim[i + j*vertNum] = m_dampingFactor * IN + 1 - m_dampingFactor;
                }
                else {
                	float OUT = minOutDeg / maxOutDeg;
                    mSim[i + j*vertNum] = m_dampingFactor*((1 - ioBalFactor)*IN + ioBalFactor * OUT) + 1 - m_dampingFactor;
                }
            } // end of if
        }
    } // end of outer for

} // end of initialise()



void DegRatioInitAlgor::initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
		float* mSim, float ioBalFactor, bool* mbCompEntries)
{
	initialise(vvInNeigh, vvOutNeigh, mSim, ioBalFactor);

	int vertNum = vvInNeigh.size();

    for (int i = 0; i < vertNum; ++i) {
    	for (int j = i + 1; j < vertNum; ++j) {
    		mbCompEntries[i + j * vertNum] = true;
    		mbCompEntries[j + i * vertNum] = true;
    	}
    }

    // initialise all diagonal elements to false
    for (int i = 0; i < vertNum; ++i) {
    	mbCompEntries[i + i * vertNum] = false;
    } // end of outer for
} // end of initialise()



void DegRatioInitAlgor::initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
		C_INT_PAIR_HASH_SET* hPairSim, C_INT_PAIR_HASH_MAP* vPrevSim, C_INT_PAIR_HASH_MAP* vCurrSim, float ioBalFactor)
{
	using namespace std;

	for (typename C_INT_PAIR_HASH_SET::iterator pit = hPairSim->begin(); pit != hPairSim->end(); ++pit) {
		std::pair<int, int> key = *pit;
		int i = key.first, j = key.second;
//		int i = pit->first;
//		int j = pit->second;

    	int inSizei = vvInNeigh[i].size();
    	int inSizej = vvInNeigh[j].size();
    	int outSizei = vvOutNeigh[i].size();
    	int outSizej = vvOutNeigh[j].size();
        int maxInDeg = max(inSizei, inSizej);
        int maxOutDeg = max(outSizei, outSizej);
        float minInDeg = min(inSizei, inSizej);
        float minOutDeg = min(outSizei, outSizej);

//        std::pair<int, int> key(i,j);
//        int currIndex = pit->second;
        if (maxInDeg == 0) {
            if (maxOutDeg == 0) {
//            	hPairSim.insert(make_pair(key, 0));
            	vPrevSim->insert(make_pair(key, 0));
            	// might not be necessary, but do it anyway to avoid future bug issues
            	vCurrSim->insert(make_pair(key, 0));
            }
            else {
                float OUT = minOutDeg / maxOutDeg;
//                hPairSim.insert(make_pair(key, m_dampingFactor * OUT + 1 -s m_dampingFactor));
                vPrevSim->insert(make_pair(key, m_dampingFactor * OUT + 1 - m_dampingFactor));
                vCurrSim->insert(make_pair(key, m_dampingFactor * OUT + 1 - m_dampingFactor));
            }
        }
        else {
            float IN = minInDeg / maxInDeg;
            if (maxOutDeg == 0) {
//            	hPairSim.insert(make_pair(key, m_dampingFactor * IN + 1 -s m_dampingFactor));
            	vPrevSim->insert(make_pair(key, m_dampingFactor * IN + 1 - m_dampingFactor));
            	vCurrSim->insert(make_pair(key, m_dampingFactor * IN + 1 - m_dampingFactor));
            }
            else {
                float OUT = minOutDeg / maxOutDeg;
//                hPairSim.insert(make_pair(key, m_dampingFactor*((1 - ioBalFactor)*IN + ioBalFactor * OUT) + 1 - m_dampingFactor));
                vPrevSim->insert(make_pair(key, m_dampingFactor*((1 - ioBalFactor)*IN + ioBalFactor * OUT) + 1 - m_dampingFactor));
                vCurrSim->insert(make_pair(key, m_dampingFactor*((1 - ioBalFactor)*IN + ioBalFactor * OUT) + 1 - m_dampingFactor));
            }
        } // end of if

	}
} // end of initialise() - 2

/* ********************************************************** */

BinaryInitAlgor::BinaryInitAlgor(float dampingFactor)
	: InitAlgor(dampingFactor)
{
}


BinaryInitAlgor::~BinaryInitAlgor()
{}

void BinaryInitAlgor::initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
		float* mSim, float ioBalFactor)
{
	using namespace std;

	int vertNum = vvInNeigh.size();

	// intialise all non-diagonal elements to 0
    for (int i = 0; i < vertNum; ++i) {
    	for (int j = i + 1; j < vertNum; ++j) {
    		mSim[i + j * vertNum] = 0;
    		mSim[j + i * vertNum] = 0;
    	}
    }

    // initialise all diagonal elements to 1
    for (int i = 0; i < vertNum; ++i) {
    	mSim[i + i * vertNum] = 1;
    } // end of outer for

} // end of initialise()


void BinaryInitAlgor::initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
		float* mSim, float ioBalFactor, bool* mbCompEntries)
{
	initialise(vvInNeigh, vvOutNeigh, mSim, ioBalFactor);

	int vertNum = vvInNeigh.size();

    for (int i = 0; i < vertNum; ++i) {
    	for (int j = i + 1; j < vertNum; ++j) {
    		mbCompEntries[i + j * vertNum] = true;
    		mbCompEntries[j + i * vertNum] = true;
    	}
    }

    // initialise all diagonal elements to false
    for (int i = 0; i < vertNum; ++i) {
    	mbCompEntries[i + i * vertNum] = false;
    } // end of outer for
} // end of initialise()




void BinaryInitAlgor::initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
		C_INT_PAIR_HASH_SET* hPairSim, C_INT_PAIR_HASH_MAP* vPrevSim, C_INT_PAIR_HASH_MAP* vCurrSim, float ioBalFactor)
{

	// since we don't store diagonal elements, we only store non-diagonals, which are all 0
	for (typename C_INT_PAIR_HASH_SET::iterator pit = hPairSim->begin(); pit != hPairSim->end(); ++pit) {
//    	std::pair<int, int> key(pit->first, pit->second);

//    	hPairSim.insert(make_pair(key, 0));
		vPrevSim->insert(std::make_pair(*pit, 0));
		vCurrSim->insert(std::make_pair(*pit, 0));
	}
} // end of initialise() - 2


/* ********************************************************** */

DegBinaryInitAlgor::DegBinaryInitAlgor(float dampingFactor)
	: InitAlgor(dampingFactor)
{

}


DegBinaryInitAlgor::~DegBinaryInitAlgor()
{}

void DegBinaryInitAlgor::initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
		float* mSim,  float ioBalFactor)
{
	using namespace std;

	int vertNum = vvInNeigh.size();

	// non-diagonal elements
    for (int i = 0; i < vertNum; ++i) {
        for (int j = i+1; j < vertNum; ++j) {
        	int inDegi = vvInNeigh[i].size();
        	int inDegj = vvInNeigh[j].size();
        	int outDegi = vvOutNeigh[i].size();
        	int outDegj = vvOutNeigh[j].size();

        	if (inDegi == inDegj && outDegi == outDegj) {
        		mSim[i + j*vertNum] = 1;
        		mSim[j + i*vertNum] = 1;
        	}
        	else {
        		mSim[i + j*vertNum] = 0;
        		mSim[j + i*vertNum] = 0;
        	}
        }
    } // end of outer for

    // diagonal elements (always 1)
    // initialise all diagonal elements to 1
    for (int i = 0; i < vertNum; ++i) {
    	mSim[i + i * vertNum] = 1;
    } // end of outer for


} // end of initialise()


void DegBinaryInitAlgor::initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
		float* mSim, float ioBalFactor, bool* mbCompEntries)
{
	initialise(vvInNeigh, vvOutNeigh, mSim, ioBalFactor);

	int vertNum = vvInNeigh.size();

    for (int i = 0; i < vertNum; ++i) {
    	for (int j = i + 1; j < vertNum; ++j) {
    		mbCompEntries[i + j * vertNum] = true;
    		mbCompEntries[j + i * vertNum] = true;
    	}
    }

    // initialise all diagonal elements to false
    for (int i = 0; i < vertNum; ++i) {
    	mbCompEntries[i + i * vertNum] = false;
    } // end of outer for
} // end of initialise()


void DegBinaryInitAlgor::initialise(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh,
		C_INT_PAIR_HASH_SET* hPairSim, C_INT_PAIR_HASH_MAP* vPrevSim, C_INT_PAIR_HASH_MAP* vCurrSim, float ioBalFactor)
{

	for (typename C_INT_PAIR_HASH_SET::iterator pit = hPairSim->begin(); pit != hPairSim->end(); ++pit) {
//		int i = pit->first;
//		int j = pit->second;
		std::pair<int, int> key = *pit;
		int i = key.first, j = key.second;

    	int inDegi = vvInNeigh[i].size();
    	int inDegj = vvInNeigh[j].size();
    	int outDegi = vvOutNeigh[i].size();
    	int outDegj = vvOutNeigh[j].size();

//    	std::pair<int, int> key(i,j);
    	int currIndex = pit->second;
    	if (inDegi == inDegj && outDegi == outDegj) {
//    		hPairSim.insert(make_pair(key, 1));
    		vPrevSim->insert(std::make_pair(key, 1));
    		vCurrSim->insert(std::make_pair(key, 1));
    	}
    	else {
//    		hPairSim.insert(make_pair(key, 0));
    		vPrevSim->insert(std::make_pair(key, 0));
    		vCurrSim->insert(std::make_pair(key, 0));
    	}
	}

} // end of initialise() - 2
