/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef PREDICTIVESUBSPACEGSGREFINEMENT_HPP_
#define PREDICTIVESUBSPACEGSGREFINEMENT_HPP_

#include "SubspaceGSGRefinement.hpp"
#include "base/grid/generation/functors/PredictiveRefinementIndicator.hpp"

namespace sg {
namespace base {

class PredictiveSubspaceGSGRefinement: public sg::base::SubspaceGSGRefinement {
public:
	PredictiveSubspaceGSGRefinement(AbstractRefinement* refinement,size_t dim): SubspaceGSGRefinement(refinement,dim){};


protected:

	void updateAdmissibleSubspaces(GridStorage* storage,
			RefinementFunctor* functor,
			ErrorVector* addedInLastRefinement,
			ErrorStorage* admissibleSubspaces);

	virtual void collectRefinableSubspaces(GridStorage* storage,
				RefinementFunctor* functor,
				HashErrorStorage* errorStorage);
};


} /* namespace base */
} /* namespace sg */
#endif /* PREDICTIVESUBSPACEGSGREFINEMENT_HPP_ */
