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
#include "base/grid/generation/refinement_strategy/RefinementDecorator.hpp"

namespace sg {
namespace base {

class PredictiveANOVARefinement: public PredictiveRefinement {
public:
	PredictiveANOVARefinement(AbstractRefinement* refinement): PredictiveRefinement(refinement){};
	virtual ~PredictiveANOVARefinement();


};

} /* namespace base */
} /* namespace sg */
#endif /* PREDICTIVEANOVAREFINEMENT_HPP_ */
