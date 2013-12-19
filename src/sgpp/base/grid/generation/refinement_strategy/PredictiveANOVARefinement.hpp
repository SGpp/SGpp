/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef PREDICTIVEANOVAREFINEMENT_HPP_
#define PREDICTIVEANOVAREFINEMENT_HPP_

#include "PredictiveRefinement.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "base/grid/generation/refinement_strategy/PredictiveRefinement.hpp";
#include "base/grid/generation/refinement_strategy/ANOVARefinement.hpp";

namespace sg {
namespace base {

class PredictiveANOVARefinement: public PredictiveRefinement, public ANOVARefinement {

public:

	PredictiveANOVARefinement(AbstractRefinement* refinement):PredictiveRefinement(refinement),ANOVARefinement(refinement), RefinementDecorator(refinement){};

	void free_refine(GridStorage* storage, RefinementFunctor* functor)
	{
		ANOVARefinement::free_refine(storage,functor);
	};
protected:

	virtual void refineGridpointsCollection(GridStorage* storage,
			RefinementFunctor* functor, size_t refinements_num, size_t* max_indices,
			RefinementFunctor::value_type* max_values)
	{
		ANOVARefinement::refineGridpointsCollection(storage,functor,refinements_num,max_indices,max_values);
	};

	virtual void collectRefinablePoints(
			GridStorage* storage, RefinementFunctor* functor,
			size_t refinements_num, size_t* max_indices,
			RefinementFunctor::value_type* max_values)
	{
		PredictiveRefinement::collectRefinablePoints(storage,functor,refinements_num,max_indices,max_values);
	};
};

} /* namespace base */
} /* namespace sg */
#endif /* PREDICTIVEANOVAREFINEMENT_HPP_ */
