/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef PREDICTIVEREFINEMENT_HPP_
#define PREDICTIVEREFINEMENT_HPP_

#include "RefinementDecorator.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "base/grid/generation/functors/PredictiveRefinementIndicator.hpp"
#include <queue>

namespace sg {
namespace base {

class PredictiveRefinement: public virtual RefinementDecorator {
public:
	PredictiveRefinement(AbstractRefinement* refinement): RefinementDecorator(refinement){};

protected:

	virtual void collectRefinablePoints(GridStorage* storage,
				RefinementFunctor* functor, size_t refinements_num, size_t* max_indices,
				RefinementFunctor::value_type* max_values);

	void free_refine(GridStorage* storage, RefinementFunctor* functor);
};

} /* namespace base */
} /* namespace sg */
#endif /* PREDICTIVEREFINEMENT_HPP_ */
