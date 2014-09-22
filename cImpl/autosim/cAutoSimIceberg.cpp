/*
 * cAutoSimIceberg.cpp
 *
 *  Created on: 03/09/2014
 *      Author: jkct440
 */



#include <iostream>
#include <utility>
#include <vector>
#include <list>
#include <unordered_set>
#include <stdexcept>
#include <cstdlib>
#include <cassert>
#include <cmath>


#include "cAutoSimIceberg.h"
#include "AutoSimInit.h"
#include "../matching/BMatching.h"
#include "../utils/matUtils.h"
#include "../utils/pairUtils.h"


#ifdef _ZERO_SIM_COUNT
	std::list<int> G_vZeroSimCount;
#endif









//AutoSimIceberg::AutoSimIceberg(float dampingFactor, int maxIter, const std::string& sInitAlgor, float simThres, float approxFaction) throw(std::invalid_argument)
//	: AutoSim(dampingFactor, maxIter, sInitAlgor), m_simThres(simThres), m_approxFaction(approxFaction)
//{
//} // end of AutoSim()
//
//
//AutoSimIceberg::AutoSimIceberg(float dampingFactor, int maxIter, float convEpsilon, const std::string& sInitAlgor, float simThres, float approxFaction) throw(std::invalid_argument)
//	: AutoSim(dampingFactor, maxIter, convEpsilon, sInitAlgor), m_simThres(simThres), m_approxFaction(approxFaction)
//{
//}



AutoSimIceberg::AutoSimIceberg(float dampingFactor, int maxIter, const std::string& sInitAlgor, bool earlySimStop, float earlySimStopThres, bool useInputBalance, float ioBalance, bool useIceberg, float simThres, float approxFaction) throw(std::invalid_argument)
	: AutoSim(dampingFactor, maxIter, sInitAlgor, earlySimStop, earlySimStopThres, useInputBalance, ioBalance), m_bIceberg(useIceberg), m_simThres(simThres), m_approxFaction(approxFaction), m_bEarlySimStop2(earlySimStop), m_earlySimStopThres2(earlySimStopThres)
{
} // end of AutoSim()


AutoSimIceberg::AutoSimIceberg(float dampingFactor, int maxIter, float convEpsilon, const std::string& sInitAlgor, bool earlySimStop, float earlySimStopThres, bool useInputBalance, float ioBalance, bool useIceberg, float simThres, float approxFaction) throw(std::invalid_argument)
	: AutoSim(dampingFactor, maxIter, convEpsilon, sInitAlgor, earlySimStop, earlySimStopThres, useInputBalance, ioBalance), m_bIceberg(useIceberg), m_simThres(simThres), m_approxFaction(approxFaction), m_bEarlySimStop2(earlySimStop), m_earlySimStopThres2(earlySimStopThres)
{
} // end of AutoSim()


AutoSimIceberg::~AutoSimIceberg()
{
} // end of ~AutoSim()





float* AutoSimIceberg::computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum) {
	using namespace std;

    assert(m_dampingFactor >= 0 && m_dampingFactor <= 1);

    // similarity matrix (column-major)
    C_INT_PAIR_HASH_MAP* pmPrevSim = new C_INT_PAIR_HASH_MAP;
    C_INT_PAIR_HASH_MAP* pmCurrSim = new C_INT_PAIR_HASH_MAP;
    // valid pair indexing
    C_INT_PAIR_HASH_SET* pmValidPairs = new C_INT_PAIR_HASH_SET;
//    C_INT_PAIR_HASH_MAP **pmPrevSim = &mPrevSim;
//    C_INT_PAIR_HASH_MAP **pmCurrSim = &mCurrSim;

    // construct neighbour list
    vector< vector<int> > vvInNeigh(vertNum);
    vector< vector<int> > vvOutNeigh(vertNum);

    // set the neighbourhoods and degrees
    std::list<int>::const_iterator sit = vSrc.begin(), tit = vTar.begin();
    for ( ; sit != vSrc.end(); ++sit, ++tit) {
//    	cout << *sit << ", " << *tit << endl;
    	assert(*sit < vertNum && *tit < vertNum);
        vvInNeigh[*tit].push_back(*sit);
        vvOutNeigh[*sit].push_back(*tit);
    } // end of for



    // parameter beta: Linear here, also can use square and Euler's
    int srcNum = 0;
    int snkNum = 0;
    for (int v = 0; v < vertNum; ++v) {
        if (vvOutNeigh[v].size() == 0 && vvInNeigh[v].size() > 0) {
            snkNum += 1;
        }
        if (vvInNeigh[v].size() == 0 && vvOutNeigh[v].size() > 0) {
            srcNum += 1;
        }
    }

    // initialise beta
    float beta = 0.5;
    if (srcNum > 0 || snkNum > 0) {
    	beta= static_cast<float>(vertNum - snkNum) / (2*vertNum - snkNum - srcNum);
    }

//    cout << "beta = " << beta << endl;
    if (m_bUseInputBalance) {
    	beta = m_ioBalance;
    }


    // compute iceberg to determine which vertex pairs
    icebergFilter(vvInNeigh, vvOutNeigh, pmValidPairs, beta);


    // initialise similarity matrix
    m_pfInitAlgor->initialise(vvInNeigh, vvOutNeigh, pmValidPairs, pmPrevSim, pmCurrSim, beta);

//    // initialise pmCurrSim to all 0's
//    for (typename C_INT_PAIR_HASH_MAP::iterator pit = pmPrevSim->begin(); pit != pmPrevSim->end(); ++pit) {
//    	pmCurrSim->insert(make_pair(pit->first, 0));
//    }


#ifdef _COLLECT_MATCHING_CHANGES_
    // resize matching matrices to right size
    m_prevInMatchingPairMatrix.resize(vertNum * vertNum);
    m_prevOutMatchingPairMatrix.resize(vertNum * vertNum);
    int changedMatchingNum = 0;
#endif


    // temporary structure for mIn and mOut
    // TODO: use max in and out degree instead
    vector<float> mIn(vertNum * vertNum);
    vector<float> mOut(vertNum * vertNum);

    // perform loop iterations
    for (int t = 1; t <= m_maxIter; ++t) {
    	cout << "iteration " << t << endl;

#ifdef _ZERO_SIM_COUNT
        int currZeroSimCount = 0;
#endif

//        C_SimVector mTempPrevSim = *pmPrevSim;
//        C_SimVector mTempCurrSim = *pmCurrSim;

        // loop through pairs
        for (typename C_INT_PAIR_HASH_SET::const_iterator currPairIt = pmValidPairs->begin(); currPairIt != pmValidPairs->end(); ++currPairIt) {

        	std::pair<int, int> key = *currPairIt;
        	int i = key.first;
        	int j = key.second;

            // In matching
            float inMatchCost = 0;
            int inDegI = vvInNeigh[i].size(), inDegJ = vvInNeigh[j].size();
            if (inDegI > 0 && inDegJ > 0) {
//                    float *mIn = new float[vInDeg[i] * vInDeg[j]];
            	// initialise the cost matrix
                for (int xi = 0; xi < inDegI; ++xi) {
                	for (int yj = 0; yj < inDegJ; ++yj) {
                		// TODO: make sure order is correct
                		std::pair<int, int> searchKey(vvInNeigh[i][xi], vvInNeigh[j][yj]);
                		typename C_INT_PAIR_HASH_MAP::const_iterator searchIt = pmPrevSim->find(searchKey);
                		if (searchIt != pmPrevSim->end()) {
                			mIn[xi + yj*inDegI] = searchIt->second;
                		}
                		// else estimate
                		else {
                			mIn[xi + yj*inDegI]	= estimateSim(vvInNeigh[i][xi], vvInNeigh[j][yj], vvInNeigh, vvOutNeigh, beta);
                		}
                    }
                }


                vector<int> m1;
                vector<int> m2;
                inMatchCost = matching(inDegI, inDegJ, mIn, m1, m2);
#ifdef _COLLECT_MATCHING_CHANGES_
                if (t > 1) {
                	bool bSame = compareMatching(m1, m2, m_prevInMatchingPairMatrix[i + j*vertNum].first, m_prevInMatchingPairMatrix[i + j*vertNum].second);
                    if (bSame) {
                    	changedMatchingNum += 1;
                    }
                }

                m_prevInMatchingPairMatrix[i + j*vertNum] = make_pair(m1, m2);
#endif
                assert(inDegI > 0 && inDegJ > 0);
                inMatchCost = inMatchCost / max(inDegI, inDegJ);
            }


            float outMatchCost = 0;
            int outDegI = vvOutNeigh[i].size(), outDegJ = vvOutNeigh[j].size();
            if (outDegI > 0 && outDegJ > 0) {
//                    float *mOut = new float[vOutDeg[i] * vOutDeg[j]];

            	// initialise the cost matrix
                for (int xi = 0; xi < outDegI; ++xi) {
                	for (int yj = 0; yj < outDegJ; ++yj) {
                		std::pair<int, int> searchKey(vvOutNeigh[i][xi], vvOutNeigh[j][yj]);
                		typename C_INT_PAIR_HASH_MAP::const_iterator searchIt = pmPrevSim->find(searchKey);
                		if (searchIt != pmPrevSim->end()) {
                			mOut[xi + yj*outDegI] = searchIt->second;
                		}
                		// else estimate
                		else {
                			mOut[xi + yj*outDegI]	= estimateSim(vvOutNeigh[i][xi], vvOutNeigh[j][yj], vvInNeigh, vvOutNeigh, beta);
                		}
                    }
                }


                vector<int> m1;
                vector<int> m2;
                outMatchCost = matching(outDegI, outDegJ, mOut, m1, m2);
#ifdef _COLLECT_MATCHING_CHANGES_
                    if (t > 1) {
                    	bool bSame = compareMatching(m1, m2, m_prevOutMatchingPairMatrix[i + j*vertNum].first, m_prevOutMatchingPairMatrix[i + j*vertNum].second);
                    	if (bSame) {
                    		changedMatchingNum += 1;
                    	}
                    }
                    m_prevOutMatchingPairMatrix[i + j*vertNum] = make_pair(m1, m2);
#endif

                assert(outDegI > 0 && outDegJ > 0);
                outMatchCost = outMatchCost / max(outDegI, outDegJ);
            }

            typename C_INT_PAIR_HASH_MAP::iterator currValueIt = pmCurrSim->find(key);
            assert(currValueIt != pmCurrSim->end());
            // compute the whole value
            if (max(inDegI, inDegJ) == 0) {
            	if (max(outDegI, outDegJ) == 0) {
            		currValueIt->second = 0;
            		//mTempCurrSim[i + j*vertNum] = 0;
                }
                else {
                	currValueIt->second = outMatchCost * m_dampingFactor + 1 - m_dampingFactor;
//                    mTempCurrSim[i + j*vertNum] = outMatchCost * m_dampingFactor + 1 - m_dampingFactor;
                }
            }
            else {
				if (max(outDegI, outDegJ) == 0) {
					currValueIt->second = inMatchCost * m_dampingFactor + 1 - m_dampingFactor;
//                        mTempCurrSim[i + j*vertNum] = inMatchCost * m_dampingFactor + 1 - m_dampingFactor;
				}
				else {
					currValueIt->second = m_dampingFactor * ((1-beta)*inMatchCost+ beta*outMatchCost) + 1 - m_dampingFactor;
//                    	mTempCurrSim[i + j*vertNum] = m_dampingFactor * ((1-beta)*inMatchCost+ beta*outMatchCost) + 1 - m_dampingFactor;
				}
            }
                // assign the other symmetric similarity
//                mTempCurrSim[j + i*vertNum] = mTempCurrSim[i + j*vertNum];

#ifdef _ZERO_SIM_COUNT
                if (mTempCurrSim[i + j*vertNum] == 0) {
                	currZeroSimCount += 2;
                }
#endif


    	} // end of loop through pairs

//        // loop through diagonal pairs
//        for (int i = 0; i < vertNum; ++i) {
//            // for all initialisation schemes used so far, S(u,u) = 1, so a maximum matching would be 1, so the
//            // similarity between i and i would also be 1
//            mTempCurrSim[i + i*vertNum] = 1;
//        } // end of outer for



#ifdef _ZERO_SIM_COUNT
             G_vZeroSimCount.push_back(currZeroSimCount);
#endif


#ifdef _COLLECT_SIM_DELTA_
    	// do comparison
        if (t > 1) {
        	assert(*pmPrevSim != NULL);
        	assert(*pmCurrSim != NULL);

        	m_vSimDelta.push_back(l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum));
        }
#endif

        // see if we analyse early stop
        if (m_bEarlySimStop && t > 1) {
        	filterEarlySimPairs(pmValidPairs, pmPrevSim, pmCurrSim);
        }


        m_iterRan = t;

        // check if we have convergence epsilon
        if (m_bUseConvEpsilon) {
//        	cout << "l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum) = " << l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum) << endl;
        	if (matDiff(pmValidPairs, pmPrevSim, pmCurrSim) < m_convEpsilon) {
        		break;
        	}
        }


        // swap the sim matrices so we do not need to allocate more matrices then need to be
        C_INT_PAIR_HASH_MAP *pmTempSim = pmPrevSim;
        pmPrevSim = pmCurrSim;
        pmCurrSim = pmTempSim;
    } // end of loop through iterations

    // destroy dynamically allocated memory
    delete pmPrevSim;
    pmPrevSim = NULL;

    // copy to output matrix output
    // TODO: add timer so we can deduct this time from iceberg running time
    float* mSim = new float[vertNum * vertNum];

	for (int i = 0; i < vertNum; ++i) {
		for (int j = i+1; j < vertNum; ++j) {
			typename C_INT_PAIR_HASH_MAP::const_iterator pit = pmCurrSim->find(make_pair(i,j));
			if (pit != pmCurrSim->end()) {
				mSim[i + j*vertNum] =  pit->second;
				mSim[j + i*vertNum] =  pit->second;
			}
			else {
				mSim[i + j*vertNum] = estimateSim(i, j, vvInNeigh, vvOutNeigh, beta);
				mSim[j + i*vertNum] = mSim[i + j*vertNum];
			}
		}
	}

	// diagonal are all ones
	for (int i = 0; i < vertNum; ++i) {
		mSim[i + i*vertNum] = 1;
	}

	delete pmCurrSim;

    return mSim;
} // end of AutoSim::computeSim()




/* ********************************************************************* */


void AutoSimIceberg::icebergFilter(const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh, C_INT_PAIR_HASH_SET* hValidPair, float beta)
{
	using namespace std;

	int vertNum = vvInNeigh.size();
	assert(vvInNeigh.size() == vvOutNeigh.size());

	float thresRatio = static_cast<float>(m_simThres - 1 + m_dampingFactor) / m_dampingFactor;

	for (int i = 0; i < vertNum; ++i) {
		for (int j = i + 1; j < vertNum; ++j) {
			// we only do iceberg filtering if need to be
			if (m_bIceberg) {
				// use float to avoid integer devision
				float inDegi = vvInNeigh[i].size();
				int inDegj = vvInNeigh[j].size();

				float inRatio;
				if (inDegi >= inDegj) {
					inRatio = inDegj / inDegi;
				}
				else {
					inRatio = inDegi / inDegj;
				}

				// use float to avoid integer devision
				float outDegi = vvOutNeigh[i].size();
				int outDegj = vvOutNeigh[j].size();

				float outRatio;
				if (outDegi >= outDegj) {
					outRatio = outDegj / outDegi;
				}
				else {
					outRatio = outDegi / outDegj;
				}


				// perform rule1 test now
				if ((1-beta) * inRatio + beta * outRatio > thresRatio) {

					// perform rule2 test
					// find maximum matching between u and v

					hValidPair->insert(make_pair(i, j));
				}
			}
			// no iceberg so insert all pairs into hValidPair
			else {
				hValidPair->insert(make_pair(i,j));
			}

		} // end of inner for loop over vertices
	} // end of outer for loop

} // end of icebergFilter()



float AutoSimIceberg::estimateSim(int vert1, int vert2, const std::vector< std::vector<int> >& vvInNeigh, const std::vector< std::vector<int> >& vvOutNeigh, float beta) const
{

	// if same vertices, we return a sim estimate of 1
	if (vert1 == vert2) {
		return 1;
	}

	// find smaller neighbourhoods
	int inDeg1 = vvInNeigh[vert1].size();
	float inDeg2 = vvInNeigh[vert2].size();
	int outDeg1 = vvOutNeigh[vert1].size();
	float outDeg2 = vvOutNeigh[vert2].size();

	float inRatio, outRatio;
	if (inDeg1 == 0 && inDeg2 == 0) {
		inRatio = 0;
	}
	else {
		if (inDeg1 <= inDeg2) {
			assert(inDeg2 > 0);
			inRatio = inDeg1/ inDeg2;
		}
		else {
			assert(inDeg1 > 0);
			inRatio = inDeg2 / inDeg1;
		}
	}

	if (outDeg1 ==0 && outDeg2 == 0) {
		outRatio = 0;
	}
	else {
		if (outDeg1 <= outDeg2) {
			assert(outDeg2 > 0);
			outRatio = outDeg1/ outDeg2;
		}
		else {
			assert(outDeg1 > 0);
			outRatio = outDeg2 / outDeg1;
		}
	}


	return m_approxFaction * m_dampingFactor * ((1-beta) * inRatio + beta * outRatio) + 1 - m_dampingFactor;
} // end of estimateSim()



void AutoSimIceberg::filterEarlySimPairs(C_INT_PAIR_HASH_SET* phValidPair, const C_INT_PAIR_HASH_MAP* const pmPrevSim, const C_INT_PAIR_HASH_MAP* const pmCurrSim)
{
	using namespace std;

	// list of pairs to remove
	list<pair<int, int> > lToRemove;

	for (typename C_INT_PAIR_HASH_SET::const_iterator pit = phValidPair->begin(); pit != phValidPair->end(); ++pit) {
		C_INT_PAIR_HASH_MAP::const_iterator prevIt = pmPrevSim->find(*pit);
		C_INT_PAIR_HASH_MAP::const_iterator currIt = pmCurrSim->find(*pit);

		if (fabs(prevIt->second - currIt->second) <= m_earlySimStopThres2) {
			lToRemove.push_back(*pit);
		}
	}

	// remove the actual pairs
	for (list<pair<int, int> >::const_iterator lit = lToRemove.begin(); lit != lToRemove.end(); ++lit) {
		phValidPair->erase(*lit);
	}
} // end of filterEarlySimPairs()




float AutoSimIceberg::matDiff(const C_INT_PAIR_HASH_SET* const phValidPair, const C_INT_PAIR_HASH_MAP* const pmPrevSim, const C_INT_PAIR_HASH_MAP* const pmCurrSim) const
{
	float diff = 0;

	for (typename C_INT_PAIR_HASH_SET::const_iterator pit = phValidPair->begin(); pit != phValidPair->end(); ++pit) {
		C_INT_PAIR_HASH_MAP::const_iterator prevIt = pmPrevSim->find(*pit);
		C_INT_PAIR_HASH_MAP::const_iterator currIt = pmCurrSim->find(*pit);

		diff += fabs(prevIt->second - currIt->second);
	}

	// normalise with the number of valid pairs
	return diff / phValidPair->size();
} // end of l1MatDiff()

