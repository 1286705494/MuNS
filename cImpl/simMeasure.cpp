
#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <cstdlib>
#include <cstring>
#include <unordered_set>
#include <boost/tokenizer.hpp>

#include "autosim/cAutoSim.h"
#include "autosim/cAutoSimIceberg.h"
#include "rolesim/cRoleSim.h"
#include "simrank/cSimRank.h"
#include "prank/cPRank.h"
#include "rege/cRege.h"
#include "matchsim/cMatchSim.h"
#include "matching/BMatching.h"
#include "timer/TimerHandler.h"

/* ******************************************************** */

void usage(char* sProgname);static
int getOptions(int argc, char* argv[]);



// Command line options
extern char* optarg;
//extern int optind, optopt;
extern int optind;

// default values for input parameters of executable
int g_iterInfo = 2; // (need to convert to int if use as maxIter)
float g_convEpsilon = 0.1;
bool g_bUseConvEpsilon = false;

float g_dampingFactor = 0.9;
std::string g_initAlgorName = "degRatioInit";
float g_ioBalance = 0.5;
bool g_bUseInputBalance = false;
float g_icebergThres = 0.8;
float g_icebergApproxFactor = 0.5;

float g_earlySimStopThres = 0.01;


bool g_bVertSubtractOne = false;



/* ************************************************************* */


int main(int argc, char *argv[])
{
    using namespace std;
    
//    for (int a = 0; a < argc; ++a) {
//    	cout << argv[a] << endl;
//    }

	// read in parameters
	int nextOptIndex = getOptions(argc, argv);

	// there should be one or more files
	if (argc - nextOptIndex != 3) {

		cerr << "Incorrect number of arguments." << endl;
		usage(argv[0]);
	}


    char* sGraphFilename = argv[nextOptIndex];
    char* sMeasure = argv[nextOptIndex+1];
//    int vertNum = atoi(argv[nextOptIndex+2]);
    char* sSimOutFilename = argv[nextOptIndex+2];
    

    // Read in file
    ifstream fIn(sGraphFilename);
    
    
    boost::char_separator<char> sSep(",\t ");
    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
            
    unordered_set<int> vUniqueVerts;

    list<int> vSrc;
    list<int> vTar;
    

    string sLine;
    while(getline(fIn, sLine)) {
    	// we ignore lines that start with '#'
    	if (sLine[0] == '#') {
    		continue;
    	}
        // tokenize string
        Tokenizer tok(sLine, sSep);
        typename Tokenizer::const_iterator it = tok.begin();
        string sSrc = *it;
        int src = atoi(sSrc.c_str());
        ++it;
        string sTar = *it;
        int tar = atoi(sTar.c_str());

        if (g_bVertSubtractOne) {
        	src -= 1;
        	tar -= 1;
        }

        vSrc.push_back(src);
        vTar.push_back(tar);

        vUniqueVerts.insert(src);
        vUniqueVerts.insert(tar);
    }
    
    fIn.close();

    // number of vertices
    int vertNum = vUniqueVerts.size();

#ifdef _TIMING_ALGOR_
    // start timer
    TimerHandler timer;
    timer.startTimer();
#endif
    
    IterSimilarity* pfSim = NULL;
    if (strcmp(sMeasure, "simrank") == 0) {
    	if (g_bUseConvEpsilon) {
    		pfSim = new SimRank(g_dampingFactor, g_iterInfo, g_convEpsilon);
    	}
    	else {
    		pfSim = new SimRank(g_dampingFactor, g_iterInfo);
    	}
    }
    else if (strcmp(sMeasure, "prank") == 0) {
    	if (g_bUseConvEpsilon) {
    		pfSim = new PRank(g_dampingFactor, g_iterInfo, g_convEpsilon, g_ioBalance);
    	}
    	else {
    		pfSim = new PRank(g_dampingFactor, g_iterInfo, g_ioBalance);
    	}
    }
    else if (strcmp(sMeasure, "matchsim") == 0) {
    	if (g_bUseConvEpsilon) {
        	pfSim = new MatchSim(g_iterInfo, g_convEpsilon);
    	}
    	else {
    		pfSim = new MatchSim(g_iterInfo);
    	}
    }
    else if (strcmp(sMeasure, "rolesim") == 0) {
    	if (g_bUseConvEpsilon) {
    		pfSim = new RoleSim(g_dampingFactor, g_iterInfo, g_convEpsilon, g_initAlgorName);
    	}
    	else {
    		pfSim = new RoleSim(g_dampingFactor, g_iterInfo, g_initAlgorName);
    	}
    }
    else if (strcmp(sMeasure, "rege") == 0) {
    	if (g_bUseConvEpsilon) {
			pfSim = new Rege(g_iterInfo, g_convEpsilon);
		}
		else {
			pfSim = new Rege(g_iterInfo);
		}
    }
    else if (strcmp(sMeasure, "autosim") == 0) {
    	// not using earlyStop
    	if (g_bUseConvEpsilon) {
    		pfSim = new AutoSim(g_dampingFactor, g_iterInfo, g_convEpsilon, g_initAlgorName, false, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance, false);
    	}
    	// not using convergence epsilon
    	else {
    		pfSim = new AutoSim(g_dampingFactor, g_iterInfo, g_initAlgorName, false, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance, false);
    	}
    }
    else if (strcmp(sMeasure, "vertBalAutosim") == 0) {
        	// not using earlyStop
        	if (g_bUseConvEpsilon) {
        		pfSim = new AutoSim(g_dampingFactor, g_iterInfo, g_convEpsilon, g_initAlgorName, false, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance, true);
        	}
        	// not using convergence epsilon
        	else {
        		pfSim = new AutoSim(g_dampingFactor, g_iterInfo, g_initAlgorName, false, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance, true);
        	}
        }
    else if (strcmp(sMeasure, "earlyStopAutosim") == 0) {
//    	if (g_bUseConvEpsilon) {
//    		pfSim = new AutoSim(g_dampingFactor, g_iterInfo, g_convEpsilon, g_initAlgorName, true, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance);
//    	}
//    	// not using convergence epsilon
//    	else {
//    		pfSim = new AutoSim(g_dampingFactor, g_iterInfo, g_initAlgorName, true, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance);
//    	}
		if (g_bUseConvEpsilon) {
			pfSim = new AutoSimIceberg(g_dampingFactor, g_iterInfo, g_convEpsilon, g_initAlgorName, true, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance, false, g_icebergThres, g_icebergApproxFactor);
		}
		else {
			pfSim = new AutoSimIceberg(g_dampingFactor, g_iterInfo, g_initAlgorName, true, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance, false, g_icebergThres, g_icebergApproxFactor);
		}
    }
    else if (strcmp(sMeasure, "icebergAutosim") == 0) {
    	// not using earlyStop
    	if (g_bUseConvEpsilon) {
    		pfSim = new AutoSimIceberg(g_dampingFactor, g_iterInfo, g_convEpsilon, g_initAlgorName, false, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance, true, g_icebergThres, g_icebergApproxFactor);
    	}
    	else {
    		pfSim = new AutoSimIceberg(g_dampingFactor, g_iterInfo, g_initAlgorName, false, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance, true, g_icebergThres, g_icebergApproxFactor);
    	}
    }
    else if (strcmp(sMeasure, "icebergEarlyStopAutosim") == 0) {
		if (g_bUseConvEpsilon) {
			pfSim = new AutoSimIceberg(g_dampingFactor, g_iterInfo, g_convEpsilon, g_initAlgorName, true, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance, true, g_icebergThres, g_icebergApproxFactor);
		}
		else {
			pfSim = new AutoSimIceberg(g_dampingFactor, g_iterInfo, g_initAlgorName, true, g_earlySimStopThres, g_bUseInputBalance, g_ioBalance, true, g_icebergThres, g_icebergApproxFactor);
		}
	}
    else {
    	cerr << "Unknown similarity measure option " << sMeasure << endl;
    	usage(argv[0]);
    }

    // construct distances
//    double* mSim = roleSim2(vSrc, vTar, vSrc.size(), vertNum, iterNum, dampingFactor);
    assert(pfSim != NULL);

    float* mSim = pfSim->computeSim(vSrc, vTar, vSrc.size(), vertNum);


    // output
    ofstream fOut(sSimOutFilename);


    // if we are collecting delta sim information
#ifdef _COLLECT_SIM_DELTA_
    const std::vector<float>& vSimDelta = pfSim->getSimDelta();
    fOut << "SimDelta: ";
    typename std::vector<float>::const_iterator dit = vSimDelta.begin();
    if (dit != vSimDelta.end()) {
    	fOut << *dit;
    	++dit;
    }
    for ( ; dit != vSimDelta.end(); ++dit) {
    	fOut << ", " << *dit;
    }
    fOut << endl;
#endif


#ifdef _COLLECT_EARLYSTOP_STATS_
    // aggregate the statistics
    std::vector<int> vIterCount(pfSim->getIterRan(), 0);
    const std::vector<int>& mIterStopped = pfSim->getIterConverged();
    for (int i = 0; i < vertNum; ++i) {
    	for (int j = i+1; j < vertNum; ++j) {
    		++vIterCount[mIterStopped[i + j*vertNum]];
    	}
    }

    for (typename std::vector<int>::const_iterator vit = vIterCount.begin(); vit != vIterCount.end(); ++vit) {
    	fOut << *vit << ", ";
    }
    fOut << endl;
#endif


    // if we are collecting total running time information
#ifdef _TIMING_ALGOR_
    // stop timer
    timer.stopTimer();
    // output to stdout

    fOut << "RunningTime: " << timer.getElapsedTime() << endl;
#endif

    fOut << "Iteration: " << pfSim->getIterRan() << endl;

    // output upper diagonal, minus the diagonal, which is always 1 for all our measures
//	for (int i = 0; i < vertNum; ++i) {
//		for (int j = i+1; j < vertNum-1; ++j) {
//			fOut << mSim[i + j*vertNum] << ",";
//		}
//		fOut << mSim[i + (vertNum-1) * vertNum] << endl;
//	}

	for (int i = 0; i < vertNum; ++i) {
		for (int j = 0; j < vertNum-1; ++j) {
			fOut << mSim[i + j*vertNum] << ",";
		}
		fOut << mSim[i + (vertNum-1) * vertNum] << endl;
	}

    fOut.close();

    delete pfSim;

    
#ifdef _ZERO_SIM_COUNT
    cout << "zero counts" << endl;
    for (typename std::list<int>::const_iterator zit = G_vZeroSimCount.begin(); zit != G_vZeroSimCount.end(); ++zit) {
    	cout << *zit << ",";
    }
    cout << endl;
#endif



    // release memory from roleSim2()
    delete[] mSim;
} // end of main()


/* *************************************************************************** */


/*
 * Get the options, and modify the locally global flags as appropriate.
 */
static
int getOptions(int argc, char* argv[])
{
	// check if there are any options
	if (argc <= 1) {
		usage(argv[0]);
	}

	// now read in some options
	int currOpt;
	/*
	* f -
	*/
	const char* optString = "t:d:i:e:b:c:a:s:m";
	while ((currOpt = getopt(argc, argv, optString)) != -1) {
		switch (currOpt) {
			case 't':
				g_iterInfo = atoi(optarg);
				break;
			case 'e':
				g_convEpsilon = atof(optarg);
				g_bUseConvEpsilon = true;
				break;
			case 'd':
				g_dampingFactor = atof(optarg);
				break;
			case 'i':
				g_initAlgorName = optarg;
				break;
			case 'b':
				g_ioBalance = atof(optarg);
				g_bUseInputBalance = true;
				break;
			case 'c':
				g_icebergThres = atof(optarg);
				assert(g_icebergThres >= 0 && g_icebergThres <= 1);
				break;
			case 'a':
				g_icebergApproxFactor = atof(optarg);
				break;
			case 's':
				g_earlySimStopThres = atof(optarg);
				break;
			case 'm':
				g_bVertSubtractOne = true;
				break;
			default:
				std::cerr << currOpt << " is not a valid option." << std::endl;
				// error
				usage(argv[0]);
				// TODO: throw some exception rather than exit in usage.
				break;

		} // end of switch
	} // end of while

	return optind;
} // end of getOptions()



void usage(char* sProgname)
{
	using namespace std;
	cerr << sProgname << ": [inFilename] [outFilename]" << endl;
	exit(1);
}
