/*
 * CombiAtanSpecialStretching.cpp
 *
 *  Created on: Apr 5, 2011
 *      Author: benk
 */

#include "CombiAtanSpecialStretching.hpp"
#include <math.h>

using namespace std;

void combigrid::AtanSpecialStretching::get1DStretching(
		int level , double min, double max,
		std::vector<double>& stretching) const {

	int nrPoints = combigrid::powerOfTwo[level] + 1;
	stretching.resize(nrPoints);
	std::vector<double> tmpPoints(nrPoints);

	for (int ii = 0 ; ii < nrPoints; ii++)
		tmpPoints[ii] = (double)(2*ii - combigrid::powerOfTwo[level]) / (double)(combigrid::powerOfTwo[level]);

	for (int ii = 0 ; ii < nrPoints; ii++)
		stretching[ii] = ::pow(tmpPoints[ii],5) + 0.1 * tmpPoints[ii];
    // calculate grading
    for (int ii = 0 ; ii < nrPoints; ii++)
    	stretching[ii] = 3.0 * ::atan( stretching[ii] );
    // do the scaling
    for (int ii = 0 ; ii < nrPoints; ii++)
    	stretching[ii] = min + (max-min) * 0.5 * ( 1 + (stretching[ii] / stretching[nrPoints-1]));
}
