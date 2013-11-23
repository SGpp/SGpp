/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef SUBSPACECOARSENING_HPP_
#define SUBSPACECOARSENING_HPP_

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/CoarseningFunctor.hpp"
#include "base/grid/generation/hashmap/dataStructures/SubspaceCoarseningErrorContainer.hpp"
#include "base/grid/generation/hashmap/HashCoarsening.hpp"
#include "base/datatypes/DataVector.hpp"


namespace sg {
namespace base {

class SubspaceCoarsening: public HashCoarsening {

public:
	void free_coarsen(GridStorage* storage,CoarseningFunctor* functor, DataVector* alpha);

protected:
	void collectSubspacesToCoarsen (GridStorage* gridStorage, CoarseningFunctor* functor, SubspaceCoarseningErrorStorage* errorStorage);
	void pickLowestContributionSubspaces(SubspaceCoarseningErrorStorage* errorStorage,
			size_t removementsNum,
			SubspaceCoarseningErrorContainer* minContribSubspaces );
	void cleanUpErrorStorage(SubspaceCoarseningErrorStorage* errorStorage);

private:

	size_t getMaxErrorElem(SubspaceCoarseningErrorContainer* minContribSubspaces, size_t removementsNum);
};

} /* namespace base */
} /* namespace sg */
#endif /* SUBSPACECOARSENING_HPP_ */
