/*
 * CombiTS_CT.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */

#ifndef COMBITS_CT_HPP_
#define COMBITS_CT_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"
#include "CombiCombiSchemeBasis.hpp"

namespace combigrid {

/** class of the two scale combination scheme (square root CT) <br>*/
class TS_CT : public CombiSchemeBasis{

public:

	/** Ctor
	 * @param dim dimension of the scheme
	 * @param level global level */
	TS_CT( int dim , int level ) : CombiSchemeBasis(dim) {
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

	/** Ctor
	 * @param dim
	 * @param levels the level vector for the dimension adaptive case */
	TS_CT( int dim , const std::vector<int>& levels ) : CombiSchemeBasis(dim) {
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

	/** Ctor for cases when in specific dimensions no combi should be done
	 * @param dim
	 * @param levels the level vector for the dimension adaptive case
	 * @param makeCombiInDimension */
	TS_CT( int dim , const std::vector<int>& levels ,
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

private:

};
}

#endif /* COMBITS_CT_HPP_ */
