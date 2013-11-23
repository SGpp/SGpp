/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef SUBSPACEERRORMAP_HPP_
#define SUBSPACEERRORMAP_HPP_

#include <map>
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "ErrorContainer.hpp"

namespace sg {
namespace base {

struct compareLevels
{

	bool operator()(AbstractRefinement::index_type gridPointA, AbstractRefinement::index_type gridPointB)
	{
		size_t dim = gridPointA.dim();
		//std::cout<< "comparing " << gridPointA.toString() << " & " << gridPointB.toString();
		//iterate over all dimensions
		for(size_t i = 0; i< dim; ++i)
		{
			//compare level in each dimension
			//std::cout << gridPointA.getLevel(i) << "<" << gridPointB.getLevel(i) << "? , ";
			if (gridPointA.getLevel(i) < gridPointB.getLevel(i)){
				//std::cout << " - SMALLER\n";
				return true;
			}else if (gridPointA.getLevel(i) > gridPointB.getLevel(i)) {
				//std::cout << " - LARGER\n";
				return false;
			}
		}
		//std::cout << " - LARGER\n";
		return false;
	}
};

typedef std::map<AbstractRefinement::index_type,ErrorContainer,compareLevels> SubspaceErrorStorage;

} /* namespace base */
} /* namespace sg */
#endif /* SUBSPACEERRORMAP_HPP_ */
