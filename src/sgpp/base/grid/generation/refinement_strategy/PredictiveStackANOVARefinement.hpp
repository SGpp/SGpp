/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef PREDICTIVESTACKANOVAREFINEMENT_HPP_
#define PREDICTIVESTACKANOVAREFINEMENT_HPP_

#include <stddef.h>
#include <vector>
#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "base/grid/generation/refinement_strategy/dataStructures/ErrorStorage.hpp"
#include "base/grid/generation/refinement_strategy/RefinementDecorator.hpp"
#include "base/grid/generation/functors/PredictiveRefinementIndicator.hpp"

namespace sg {
namespace base {

class PredictiveStackANOVARefinement: public RefinementDecorator {
public:
	typedef std::vector<index_type> GridPointVector;
	typedef std::vector<ErrorType> ErrorVector;


	PredictiveStackANOVARefinement(AbstractRefinement* refinement, size_t dim):RefinementDecorator(refinement), availableGridPoints(dim), firstRefinement(true){};

	void free_refine(GridStorage* storage, RefinementFunctor* functor);

protected:

	void updateAdmissiblePoints(GridStorage* storage,
			RefinementFunctor* functor,
			ErrorVector* addedInLastRefinement,
			ErrorStorage* admissibleSubspaces);

	virtual void collectRefinablePoints(GridStorage* storage,
				RefinementFunctor* functor,
				HashErrorStorage* errorStorage);

	virtual void refineGridpointsCollection(GridStorage* storage,
			ErrorStorage* errorStorage,
			ErrorVector* addedInLastStep,
			RefinementFunctor* functor);

	void resetParentLeafs(GridStorage* storage, index_type* index);

	ErrorStorage availableGridPoints;
	ErrorVector addedInLastRefinement;
	bool firstRefinement;


};

} /* namespace base */
} /* namespace sg */
#endif /* PREDICTIVESTACKANOVAREFINEMENT_HPP_ */
