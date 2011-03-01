/*
 * CombiTS_CT.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */

#ifndef COMBITS_CT_HPP_
#define COMBITS_CT_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"
#include "CombiSchemeBasis.hpp"

namespace combigrid {

/** class of the two scale combination scheme (square root CT) <br>*/
class TS_CT : public CombiSchemeBasis{

public:

	/** Ctor
	 * @param dim dimension of the scheme
	 * @param level global level */
	TS_CT( int dim , int level );

	/** Ctor
	 * @param dim
	 * @param levels the level vector for the dimension adaptive case */
	TS_CT( int dim , const std::vector<int>& levels );

	/** Ctor for cases when in specific dimensions no combi should be done
	 * @param dim
	 * @param levels the level vector for the dimension adaptive case
	 * @param makeCombiInDimension */
	TS_CT( int dim , const std::vector<int>& levels ,
	 const std::vector<bool>& makeCombiInDimension );


	/** Ctor for manual steared TS scheme where the user specifies the higher and the lower levels
	 * @param dim
	 * @param minlevels the min levels
	 * @param maxlevels the max levels  */
	TS_CT( int dim , const std::vector<int>& minlevels ,
			 const std::vector<int>& maxlevels  );

private:

};
}

#endif /* COMBITS_CT_HPP_ */
