/*
 * CombiUniformStretching.cpp
 *
 *  Created on: Apr 4, 2011
 *      Author: benk
 */

#include "CombiUniformStretching.hpp"


void combigrid::UniformStretching::get1DStretching(
		int level , double min, double max,
		std::vector<double>& stretching) const{

	int nrPoints = combigrid::powerOfTwo[level]+1;
	stretching.resize( nrPoints , 0.0 );

	// set each point just uniformly
	for (int i = 0 ; i < nrPoints ; i++){
		stretching[i] = min + ((double)i)*(max-min)/(nrPoints-1);
	}
}
