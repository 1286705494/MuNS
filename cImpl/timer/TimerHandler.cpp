//
// C++ Implementation: TimerHandler.cpp
//
// Description: Implementation of TimeHandler.h
//
//
// Author: Jeffrey Chan <jkcchan@csse.unimelb.edu.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <assert.h>
#include "TimerHandler.h"

using namespace std;


TimerHandler::TimerHandler() {
} // end of TimerHandler()



void TimerHandler::startTimer() {
	// != 0 if failred
	if (gettimeofday(&startTv, NULL) != 0) {
		// TODO: Throw Exception!
		cerr << "Timer failed." << endl;
		exit(EXIT_FAILURE);
	}	
} // end of startTimer()


void TimerHandler::stopTimer() {
	if (gettimeofday(&endTv, NULL)) {
		// TODO: Throw Exception!
		cerr << "Ending Timer failed." << endl;
		exit(EXIT_FAILURE);
	}
} // end of endTimer()


double TimerHandler::getElapsedTime() const {
	// else, calculate timing
	time_t diffSec = endTv.tv_sec - startTv.tv_sec;
	suseconds_t diffMsec = 0;
	if (endTv.tv_usec < startTv.tv_usec) {
		assert(diffSec > 0);
		diffMsec = endTv.tv_usec - startTv.tv_usec + 1000000;
		diffSec -= 1;
	}
	else 
		diffMsec = endTv.tv_usec - startTv.tv_usec;
				
	return diffSec + static_cast<double>(diffMsec) / 1000000; 	
} // end of getElapsedTime()
