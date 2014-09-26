
#include <iostream>
#include <vector>
#include <list>
#include <stdexcept>
#include <cstdlib>
#include <cassert>

#include "cAutoSim.h"
#include "AutoSimInit.h"
#include "../matching/BMatching.h"
#include "../utils/matUtils.h"


#ifdef _ZERO_SIM_COUNT
	std::list<int> G_vZeroSimCount;
#endif


inline int index(int r, int c, int vertNum) { return r < c ? r + c * vertNum : c + r * vertNum; }




//AutoSim::AutoSim(float dampingFactor, int maxIter, const std::string& sInitAlgor) throw(std::invalid_argument)
//	: DampedSimilarity(dampingFactor, maxIter), m_pfInitAlgor(NULL), m_bUseInputBalance(false), m_ioBalance(0.5)
//{
//	if (sInitAlgor == "binaryInit") {
//		m_pfInitAlgor = new BinaryInitAlgor(dampingFactor);
//	}
//	else if (sInitAlgor == "degBinaryInit") {
//		m_pfInitAlgor = new DegBinaryInitAlgor(dampingFactor);
//	}
//	else if (sInitAlgor == "degRatioInit") {
//		m_pfInitAlgor = new DegRatioInitAlgor(dampingFactor);
//	}
//	else {
//		throw std::invalid_argument("AutoSim: invalid initialisation algorithm");
//	}
//} // end of AutoSim()
//
//
//AutoSim::AutoSim(float dampingFactor, int maxIter, float convEpsilon, const std::string& sInitAlgor) throw(std::invalid_argument)
//	: DampedSimilarity(dampingFactor, maxIter, convEpsilon), m_pfInitAlgor(NULL), m_bUseInputBalance(false), m_ioBalance(0.5)
//{
//	if (sInitAlgor == "binaryInit") {
//		m_pfInitAlgor = new BinaryInitAlgor(dampingFactor);
//	}
//	else if (sInitAlgor == "degBinaryInit") {
//		m_pfInitAlgor = new DegBinaryInitAlgor(dampingFactor);
//	}
//	else if (sInitAlgor == "degRatioInit") {
//		m_pfInitAlgor = new DegRatioInitAlgor(dampingFactor);
//	}
//	else {
//		throw std::invalid_argument("AutoSim: invalid initialisation algorithm");
//	}
//}



AutoSim::AutoSim(float dampingFactor, int maxIter, const std::string& sInitAlgor, bool earlySimStop, float earlySimStopThres, bool useInputBalance, float ioBalance) throw(std::invalid_argument)
	: DampedSimilarity(dampingFactor, maxIter), m_pfInitAlgor(NULL), m_bUseInputBalance(useInputBalance), m_ioBalance(ioBalance), m_bEarlySimStop(earlySimStop), m_earlySimStopThres(earlySimStopThres)
{
//	// override default constructor
//	m_bUseInputBalance = true;
//	m_ioBalance = ioBalance;
	if (sInitAlgor == "binaryInit") {
		m_pfInitAlgor = new BinaryInitAlgor(dampingFactor);
	}
	else if (sInitAlgor == "degBinaryInit") {
		m_pfInitAlgor = new DegBinaryInitAlgor(dampingFactor);
	}
	else if (sInitAlgor == "degRatioInit") {
		m_pfInitAlgor = new DegRatioInitAlgor(dampingFactor);
	}
	else {
		throw std::invalid_argument("AutoSim: invalid initialisation algorithm");
	}
} // end of AutoSim()


AutoSim::AutoSim(float dampingFactor, int maxIter, float convEpsilon, const std::string& sInitAlgor, bool earlySimStop, float earlySimStopThres,  bool useInputBalance, float ioBalance) throw(std::invalid_argument)
	: DampedSimilarity(dampingFactor, maxIter, convEpsilon), m_pfInitAlgor(NULL), m_bUseInputBalance(useInputBalance), m_ioBalance(ioBalance), m_bEarlySimStop(earlySimStop), m_earlySimStopThres(earlySimStopThres)
{
//	// override default constructor
//	m_bUseInputBalance = true;
//	m_ioBalance = ioBalance;
	if (sInitAlgor == "binaryInit") {
		m_pfInitAlgor = new BinaryInitAlgor(dampingFactor);
	}
	else if (sInitAlgor == "degBinaryInit") {
		m_pfInitAlgor = new DegBinaryInitAlgor(dampingFactor);
	}
	else if (sInitAlgor == "degRatioInit") {
		m_pfInitAlgor = new DegRatioInitAlgor(dampingFactor);
	}
	else {
		throw std::invalid_argument("AutoSim: invalid initialisation algorithm");
	}
} // end of AutoSim()


AutoSim::~AutoSim()
{
	delete m_pfInitAlgor;

#ifdef _COLLECT_MATCHING_CHANGES_
//	for (typename C_MatchingMatrix::iterator iit =  m_prevInMatchingPairMatrix.begin(); iit != m_prevInMatchingPairMatrix.end(); ++iit) {
//		if (iit->first != NULL) {
//			delete[] iit->first;
//		}
//		if (iit->second != NULL) {
//			delete[] iit->second;
//		}
//	}
//
//	for (typename C_MatchingMatrix::iterator oit =  m_prevOutMatchingPairMatrix.begin(); oit != m_prevOutMatchingPairMatrix.end(); ++oit) {
//		if (oit->first != NULL) {
//			delete[] oit->first;
//		}
//		if (oit->second != NULL) {
//			delete[] oit->second;
//		}
//	}
#endif
} // end of ~AutoSim()









float* AutoSim::computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum) {
	using namespace std;

    assert(m_dampingFactor >= 0 && m_dampingFactor <= 1);

    // similarity matrix (column-major)
    float* mPrevSim = new float[vertNum*vertNum];
    float* mCurrSim = new float[vertNum*vertNum];
    float **pmPrevSim = &mPrevSim;
    float **pmCurrSim = &mCurrSim;

//    double *mTempPrevSim1 = *pmPrevSim;
//    for (int i = 0; i < vertNum*vertNum; ++i) {
//        mTempPrevSim1[i] = 0;
//    }

    bool* mbCompEntry = NULL;
    if (m_bEarlySimStop) {
    	mbCompEntry = new bool[vertNum*vertNum];
    }




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



    // parameter b: Linear here, also can use square and Euler's
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

    // initialise b
    double b = 0.5;
    if (srcNum > 0 || snkNum > 0) {
    	b= static_cast<double>(vertNum - snkNum) / (2*vertNum - snkNum - srcNum);
    }

//    cout << "b = " << b << endl;
    if (m_bUseInputBalance) {
    	b = m_ioBalance;
    }

    // initialise similarity matrix
    if (m_bEarlySimStop) {
    	m_pfInitAlgor->initialise(vvInNeigh, vvOutNeigh, *pmPrevSim, b, mbCompEntry);
    }
    else {
    	m_pfInitAlgor->initialise(vvInNeigh, vvOutNeigh, *pmPrevSim, b);
    }


#ifdef _COLLECT_MATCHING_CHANGES_
    // resize matching matrices to right size
    m_prevInMatchingPairMatrix.resize(vertNum * vertNum);
    m_prevOutMatchingPairMatrix.resize(vertNum * vertNum);
    int changedMatchingNum = 0;
#endif


    // temporary structure for mIn and mOut
    vector<float> mIn(vertNum * vertNum);
    vector<float> mOut(vertNum * vertNum);

    // perform loop iterations
    for (int t = 1; t <= m_maxIter; ++t) {
    	cout << "iteration " << t << endl;

#ifdef _ZERO_SIM_COUNT
        int currZeroSimCount = 0;
#endif

        float* mTempPrevSim = *pmPrevSim;
        float *mTempCurrSim = *pmCurrSim;
        // loop through pairs
        for (int i = 0; i < vertNum; ++i) {
            for (int j = i+1; j < vertNum; ++j) {

                if (m_bEarlySimStop) {
					// if we don't need to compute, just go to next entry
					if (!mbCompEntry[i + j*vertNum]) {
						continue;
					}
                }

                // In matching
            	float inMatchCost = 0;
                int inDegI = vvInNeigh[i].size(), inDegJ = vvInNeigh[j].size();
                if (inDegI > 0 && inDegJ > 0) {
//                    double *mIn = new double[vInDeg[i] * vInDeg[j]];

                    // initialise the cost matrix
                    for (int xi = 0; xi < inDegI; ++xi) {
                        for (int yj = 0; yj < inDegJ; ++yj) {
                            mIn[xi + yj*inDegI] = mTempPrevSim[vvInNeigh[i][xi] + vvInNeigh[j][yj] * vertNum];
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

//                    delete[] mIn;
                }


                float outMatchCost = 0;
                int outDegI = vvOutNeigh[i].size(), outDegJ = vvOutNeigh[j].size();
                if (outDegI > 0 && outDegJ > 0) {
//                    double *mOut = new double[vOutDeg[i] * vOutDeg[j]];

                    // initialise the cost matrix
                    for (int xi = 0; xi < outDegI; ++xi) {
                        for (int yj = 0; yj < outDegJ; ++yj) {
                            mOut[xi + yj*outDegI] = mTempPrevSim[vvOutNeigh[i][xi] + vvOutNeigh[j][yj] * vertNum];
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

//                    delete[] mOut;
                }

                if (max(inDegI, inDegJ) == 0) {
                    if (max(outDegI, outDegJ) == 0) {
                        mTempCurrSim[i + j*vertNum] = 0;
                    }
                    else {
                        mTempCurrSim[i + j*vertNum] = outMatchCost * m_dampingFactor + 1 - m_dampingFactor;
                    }
                }
                else {
                    if (max(outDegI, outDegJ) == 0) {
                        mTempCurrSim[i + j*vertNum] = inMatchCost * m_dampingFactor + 1 - m_dampingFactor;
                    }
                    else {
                    	mTempCurrSim[i + j*vertNum] = m_dampingFactor * ((1-b)*inMatchCost+ b*outMatchCost) + 1 - m_dampingFactor;
                    }
                }
                // assign the other symmetric similarity
                mTempCurrSim[j + i*vertNum] = mTempCurrSim[i + j*vertNum];

#ifdef _ZERO_SIM_COUNT
                if (mTempCurrSim[i + j*vertNum] == 0) {
                	currZeroSimCount += 2;
                }
#endif



            } // end of inner for
        } // end of outer for

        // loop through diagonal pairs
        for (int i = 0; i < vertNum; ++i) {
            // for all initialisation schemes used so far, S(u,u) = 1, so a maximum matching would be 1, so the
            // similarity between i and i would also be 1
            mTempCurrSim[i + i*vertNum] = 1;
        } // end of outer for



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



        if (m_bEarlySimStop && t > 1) {
        	matChange(*pmPrevSim, *pmCurrSim, vertNum, vertNum, mbCompEntry, m_earlySimStopThres);
        }


        m_iterRan = t;

        // check if we have convergence epsilon
        if (m_bUseConvEpsilon) {
//        	cout << "l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum) = " << l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum) << endl;
        	if (l1MatDiff(*pmPrevSim, *pmCurrSim, vertNum, vertNum) / (vertNum * vertNum) < m_convEpsilon) {
        		break;
        	}
        }



        // swap the sim matrices so we do not need to allocate more matrices then need to be
        float **pmTempSim = pmPrevSim;
        pmPrevSim = pmCurrSim;
        pmCurrSim = pmTempSim;
    } // end of loop through iterations

    // destroy dynamically allocated memory
    delete[] *pmPrevSim;
    pmPrevSim = NULL;

    if (m_bEarlySimStop && mbCompEntry != NULL) {
    	delete[] mbCompEntry;
    	mbCompEntry = NULL;
    }


    return *pmCurrSim;
} // end of AutoSim::computeSim()



bool AutoSim::compareMatching(const std::vector<int>& currM1, const std::vector<int>& currM2, const std::vector<int>& prevM1, const std::vector<int>& prevM2)
{
	bool bSame = false;



	return bSame;
} // end of compareMatching





    

  






