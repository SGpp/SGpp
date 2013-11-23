/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef SUBSPACECOARSENINGERRORCONTAINER_HPP_
#define SUBSPACECOARSENINGERRORCONTAINER_HPP_

#include <map>
#include "base/grid/generation/refinement_strategy/dataStructures/SubspaceErrorMap.hpp"

namespace sg {
namespace base {

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
