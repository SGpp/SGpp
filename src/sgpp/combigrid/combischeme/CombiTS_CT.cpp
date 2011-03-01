/*
 * CombiTS_CT.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */

#include "combigrid/combischeme/CombiTS_CT.hpp"

using namespace combigrid;
using namespace std;



TS_CT::TS_CT( int dim , int level ) : CombiSchemeBasis(dim) {
	std::vector<int> level_small( dim , level/2 );
	std::vector<int> level_tmp;
	// add all the spaces for all dimensions
	for (int i = 0; i < dim ; i++){
		level_tmp = level_small;
		level_tmp[i] = level;
		levels_vector_.push_back(level_tmp);
		cofficients_.push_back(1.0);
	}
	// add the smallest space
	levels_vector_.push_back(level_small);
	cofficients_.push_back(1.0 - dim);
}



TS_CT::TS_CT( int dim , const std::vector<int>& levels ) : CombiSchemeBasis(dim) {
	std::vector<int> level_small( dim , 0 );
	std::vector<int> level_tmp;
	// calculate
	for (int i = 0; i < dim ; i++){
		level_small[i] = levels[i] / 2;
	}
	// add all the spaces for all dimensions
	for (int i = 0; i < dim ; i++){
		level_tmp = level_small;
		level_tmp[i] = levels[i];
		levels_vector_.push_back(level_tmp);
		cofficients_.push_back(1.0);
	}
	// add the smallest space
	levels_vector_.push_back(level_small);
	cofficients_.push_back(1.0 - dim);
}



TS_CT::TS_CT( int dim , const std::vector<int>& levels ,
	 const std::vector<bool>& makeCombiInDimension ) : CombiSchemeBasis(dim){
	std::vector<int> level_small( dim , 0 );
	std::vector<bool> dimActive( dim , true );
	std::vector<int> level_tmp;
	int nrActiveDimensions = dim;
	// calculate
	for (int i = 0; i < dim ; i++){
		level_small[i] = levels[i] / 2;
		if (makeCombiInDimension[i] == false ) {
			dimActive[i] = false;
			nrActiveDimensions--;
		}
	}
	// add all the spaces for all dimensions
	for (int i = 0; i < dim ; i++){
		if (dimActive[i]){
			level_tmp = level_small;
			level_tmp[i] = levels[i];
			levels_vector_.push_back(level_tmp);
			cofficients_.push_back( 1.0 );
		}
	}
	// add the smallest space
	levels_vector_.push_back(level_small);
	cofficients_.push_back(1.0 - nrActiveDimensions);
}



TS_CT::TS_CT( int dim , const std::vector<int>& minlevels ,
		 const std::vector<int>& maxlevels  ) : CombiSchemeBasis(dim){

	std::vector<int> level_small( dim , 0 );
	std::vector<int> level_tmp;
	int nrActiveDimensions = dim;

	// first create the small levels into a vector
	for (int i = 0; i < dim ; i++){
		level_small[i] = minlevels[i];
		if (minlevels[i] == maxlevels[i] ) {
			nrActiveDimensions--;
		}
	}

	// add all the spaces for all dimensions
	for (int i = 0; i < dim ; i++){
		if (minlevels[i] == maxlevels[i] ) {
			level_tmp = level_small;
			level_tmp[i] = maxlevels[i];
			levels_vector_.push_back(level_tmp);
			cofficients_.push_back( 1.0 );
		}
	}
	// add the smallest space
	levels_vector_.push_back(level_small);
	cofficients_.push_back(1.0 - nrActiveDimensions);
}
