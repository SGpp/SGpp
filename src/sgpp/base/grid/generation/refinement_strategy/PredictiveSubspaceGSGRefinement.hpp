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
	PredictiveSubspaceGSGRefinement(RefinementDecorator* decorator);

	void freeRefineSubspace(GridStorage* storage,PredictiveRefinementIndicator* errorIndicator);


protected:

	void updateAdmissibleSubspaces(GridStorage* storage,
			PredictiveRefinementIndicator* errorIndicator,
			SubspaceVector* addedInLastRefinement,
			SubspaceErrorStorage* admissibleSubspaces);

	virtual void collectRefinableSubspaces(GridStorage* storage,
				PredictiveRefinementIndicator* errorIndicator,
				SubspaceErrorStorage* errorStorage);


	virtual void refineSubspaceCollection(GridStorage* storage,
				SubspaceErrorStorage* errorStorage,
				SubspaceVector* addedInLastStep,
				PredictiveRefinementIndicator* errorIndicator,
				size_t refinements_num,
				index_type* maxErrorSubspaces,
				RefinementFunctor::value_type* maxErrorValues);
};


} /* namespace base */
} /* namespace sg */
#endif /* PREDICTIVESUBSPACEGSGREFINEMENT_HPP_ */
