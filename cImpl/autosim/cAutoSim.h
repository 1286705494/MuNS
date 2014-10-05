#include <list>
#include <string>
#include <stdexcept>

#include "../similarity/dampedSimilarity.h"
#include "AutoSimInit.h"

#ifndef _CROLE_SIM2_H_
#define _CROLE_SIM2_H_
//
//extern
//double* roleSim2(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum, int iteration,
//		double dampingFactor);


#ifdef _ZERO_SIM_COUNT
	extern
	std::list<int> G_vZeroSimCount;
#endif


class AutoSim : public DampedSimilarity
{
protected:



	/** Initialisation algorithm. */
	InitAlgor* m_pfInitAlgor;

	/** Use this beta value (in/out balance) instead. */
	bool m_bUseInputBalance;
	float m_ioBalance;

	/** Whether to use individual level balance. Note if this is true, m_bUseInputBalance is not considered. */
	bool m_bCompIndividualBalance;

	/** Use early similarity stop. */
	bool m_bEarlySimStop;
	float m_earlySimStopThres;

protected:

#ifdef _COLLECT_MATCHING_CHANGES_
	typedef std::vector<std::pair<std::vector<int>, std::vector<int> > > C_MatchingMatrix;
	C_MatchingMatrix m_prevInMatchingPairMatrix;
	C_MatchingMatrix m_prevOutMatchingPairMatrix;
	std::vector<int> m_matchingChanges;
#endif




public:

//	AutoSim(float dampingFactor, int maxIter, const std::string& sInitAlgor) throw(std::invalid_argument);
//	AutoSim(float dampingFactor, int maxIter, float convEpsilon, const std::string& sInitAlgor) throw(std::invalid_argument);

	AutoSim(float dampingFactor, int maxIter, const std::string& sInitAlgor, bool earlySimStop = false, float earlySimStopThres = 0.1, bool useInputBalance = false, float ioBalance = 0.5, float m_bCompIndividualBalance = false) throw(std::invalid_argument);
	AutoSim(float dampingFactor, int maxIter, float convEpsilon, const std::string& sInitAlgor, bool earlySimStop = false, float earlySimStopThres = 0.1, bool useInputBalance = false, float ioBalance = 0.5, float m_bCompIndividualBalance = false) throw(std::invalid_argument);

	virtual ~AutoSim();


	virtual float* computeSim(const std::list<int>& vSrc, const std::list<int>& vTar, int edgeNum, int vertNum);

protected:

	bool compareMatching(const std::vector<int>& currM1, const std::vector<int>& currM2, const std::vector<int>& prevM1, const std::vector<int>& prevM2);



}; // end of class AutoSim




#endif
