/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef SUBSPACECOARSENINGERRORCONTAINER_HPP_
#define SUBSPACECOARSENINGERRORCONTAINER_HPP_

#include <map>
#include "base/grid/generation/hashmap/HashCoarsening.hpp"

namespace sg {
namespace base {


/**
 * compareLevels used to compare keys in SubspaceCoarseningErrorStorage.
 * Two subspaces are equal if their level vectors are equal.
 */
struct compareLevels
{

	bool operator()(HashCoarsening::index_type gridPointA, HashCoarsening::index_type gridPointB)
	{
		size_t dim = gridPointA.dim();
		//iterate over all dimensions
		for(size_t i = 0; i< dim; ++i)
		{
			//compare level in each dimension
			if (gridPointA.getLevel(i) < gridPointB.getLevel(i)){
				return true;
			}else if (gridPointA.getLevel(i) > gridPointB.getLevel(i)) {
				return false;
			}
		}
		return false;
	}
};


/**
 * SubspaceCoarseningErrorContainer stores information needed for coarsening subspaces.
 */
class SubspaceCoarseningErrorContainer {
public:

	/**
	 * Constructor
	 *
	 * initicalizes error to 0 and does not allow the subspace to be coarsened.
	 */
	SubspaceCoarseningErrorContainer()
	{
		isCoarsenable = false;
		error = 0;
	}

	//gridPoints vector of the seqence numbers of all grid points from the subspace.
	std::vector<size_t> gridPoints;
	//isCoarsenable tells if subspace is admissible for removement (i.e. if all it's grid points are leaves)
	bool isCoarsenable;
	//error error indiacator - accumulated surplusses of all grid points from the subspace.
	CoarseningFunctor::value_type error;
};

// key value pair storage for coarsening subspaces. for each subspace in the grid backward's neigbourhood (key)
//there is an SubspaceCoarseningErrorContainer containing all information on admissibility and error indicator.
typedef std::map<GridStorage::index_type,SubspaceCoarseningErrorContainer,compareLevels> SubspaceCoarseningErrorStorage;

} /* namespace base */
} /* namespace sg */
#endif /* SUBSPACECOARSENINGERRORCONTAINER_HPP_ */
