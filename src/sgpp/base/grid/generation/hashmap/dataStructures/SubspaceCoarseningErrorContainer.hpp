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

struct compareLevels
{

	bool operator()(HashCoarsening::index_type gridPointA, HashCoarsening::index_type gridPointB)
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

class SubspaceCoarseningErrorContainer {
public:
	SubspaceCoarseningErrorContainer()
	{
		// TODO Auto-generated constructor stub
		isCoarsenable = false;
		error = 0;
	}

	std::vector<size_t> gridPoints;
	bool isCoarsenable;
	CoarseningFunctor::value_type error;
};

typedef std::map<GridStorage::index_type,SubspaceCoarseningErrorContainer,compareLevels> SubspaceCoarseningErrorStorage;

} /* namespace base */
} /* namespace sg */
#endif /* SUBSPACECOARSENINGERRORCONTAINER_HPP_ */
