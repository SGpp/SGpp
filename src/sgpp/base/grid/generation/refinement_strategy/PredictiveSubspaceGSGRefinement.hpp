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



/**
 * Dimension adaptive refinement, based on the GSG algorithm by Gerstner and Griebel,
 * adapted for the predictive refinement indicator.
 * Selects hierarchical subspaces where all parent subspaces exist, calculates a refinement indicator
 * for them and creates the subspaces with high error indicators. Refinable subspaces that are not created
 * are stored for later reuse.
 */
class PredictiveSubspaceGSGRefinement: public sg::base::SubspaceGSGRefinement {
public:

	/*
	 * Constructor.
	 * -Initializes an empty error storage
	 * -Sets the refinement step to be the first one.
	 *
	 * @param refinement object implementing the core functionality (e.g.
	 * refinement with or without boundaries).
	 * @param dim dimension of the grid, that should be refined.
	 */
	PredictiveSubspaceGSGRefinement(AbstractRefinement* refinement,size_t dim): SubspaceGSGRefinement(refinement,dim){};


protected:

	/**
	 * Examines the child subspaces of the subspaces from addedInLastRefinement and adds
	 * them to the error storage, if they are admissible (all grid points in the subspace have all parents in all dimensions).
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
	 * @param addedInLastRefinement a vector to store the subspaces created in this refinement step.
	 * @param errorStorage storage container for refinable subspaces, sorted by indicator value
	 */
	void updateAdmissibleSubspaces(GridStorage* storage,
			RefinementFunctor* functor,
			ErrorVector* addedInLastRefinement,
			ErrorStorage* admissibleSubspaces);


	/**
	 * Refines a grid by adding additional Subspaces according to a RefinementFunctor provided.
	 * Refines up to RefinementFunctor::getRefinementsNum() grid points if
	 * possible, and if their refinement value is larger than RefinementFunctor::start()
	 * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
	 */
	virtual void collectRefinableSubspaces(GridStorage* storage,
				RefinementFunctor* functor,
				HashErrorStorage* errorStorage);
};


} /* namespace base */
} /* namespace sg */
#endif /* PREDICTIVESUBSPACEGSGREFINEMENT_HPP_ */
