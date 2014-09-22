//
// C++ Interface: TimerHandler
//
// Description: Performs program execution timing.
//
//
// Author: Jeffrey Chan <jkcchan@csse.unimelb.edu.au>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//


#ifndef _TIMERHANDLER_H_
#define _TIMERHANDLER_H_

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>



/**
 * Handler to perform execution timing.
 */
class TimerHandler
{
private:
	// timing variables
	struct timeval startTv, endTv;
	struct timezone startTz, endTz;
	
	

public:

	/**
	 *  @param os Output stream.
	 */
	TimerHandler();


	/**
	 * Start timer.
	 */
	void startTimer();
	
	
	/**
	 * Stop timer.
	 */
	void stopTimer();
	
	
	/**
	 * Calculate time.  Assumes timer has been started and stopped.  
	 * Else, behaviour is not guaranteed.
	 */
	double getElapsedTime() const;
	
	
	/**
	 * Output to destinated output stream.
	 * @param os Output stream.
	 */
	void output(std::ostream& os = std::cout);
	
	
}; // end of class TimerHandler
	



#endif /*_TIMERHANDLER_H_*/
