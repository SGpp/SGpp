/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef SUBSPACEGSGREFINEMENT_HPP_
#define SUBSPACEGSGREFINEMENT_HPP_

#include "SubspaceRefinement.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "base/grid/generation/refinement_strategy/dataStructures/ErrorStorage.hpp"

namespace sg {
namespace base {

class SubspaceGSGRefinement: public SubspaceRefinement {
public:


	SubspaceGSGRefinement(AbstractRefinement* refinement, size_t dim);

	void freeRefineSubspace(GridStorage* storage,RefinementFunctor* functor);

protected:


	ErrorStorage availableSubspaces;
	ErrorVector addedInLastRefinement;
	bool firstRefinement;

	bool checkAdmissible(GridStorage* storage, index_type& subspace);

	virtual void refineSubspaceCollection(GridStorage* storage,
			ErrorStorage* errorStorage,
			ErrorVector* addedInLastStep,
			RefinementFunctor* functor);

	virtual void updateAdmissibleSubspaces(GridStorage* storage,
			RefinementFunctor* functor,
			ErrorVector* addedInLastRefinement,
			ErrorStorage* admissibleSubspaces);

};
} /* namespace base */
} /* namespace sg */
#endif /* SUBSPACEGSGREFINEMENT_HPP_ */
