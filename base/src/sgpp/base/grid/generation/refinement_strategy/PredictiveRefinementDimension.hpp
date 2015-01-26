// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef PREDICTIVEREFINEMENTDIMENSION_HPP_
#define PREDICTIVEREFINEMENTDIMENSION_HPP_

#include "RefinementDecorator.hpp"
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementDimensionIndicator.hpp>
#include <vector>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {


/*
 * PredictiveRefinement performs local adaptive refinement of a sparse grid using the PredictiveRefinementIndicator.
 * This way, instead of surpluses that are most often used for refinement, refinement decisions are based upon an estimation
 * to the contribution of the MSE, which is especially helpfull for regression.
 *
 */
class PredictiveRefinementDimension: public virtual RefinementDecorator {

public:
  PredictiveRefinementDimension(AbstractRefinement* refinement): RefinementDecorator(refinement){};
  void free_refine(GridStorage* storage, PredictiveRefinementDimensionIndicator* functor);

protected:

	/**
	 * Examines the grid points and stores those indices that can be refined
	 * and have maximal indicator values.
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
	 * @param refinements_num number of points to refine
	 * @param max_indices the array where the point indices should be stored
	 * @param max_values the array where the corresponding indicator values
	 * should be stored
	 */
	virtual void collectRefinablePoints(GridStorage* storage,
				RefinementFunctor* functor, size_t refinements_num, size_t* max_indices,
				PredictiveRefinementDimensionIndicator::value_type* max_values);

	virtual void refineGridpointsCollection(GridStorage* storage,
	    RefinementFunctor* functor, size_t refinements_num, size_t* max_indices,
	    PredictiveRefinementDimensionIndicator::value_type* max_values);

	virtual size_t getIndexOfMin(PredictiveRefinementDimensionIndicator::value_type* array,
	    size_t length);

};

} /* namespace base */
} /* namespace SGPP */
#endif /* PREDICTIVEREFINEMENTDIMENSION_HPP_ */