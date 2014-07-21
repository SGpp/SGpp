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




/**
 * Dimension adaptive refinement, based on the GSG algorithm by Gerstner and Griebel.
 * Selects hierarchical subspaces where all parent subspaces exist, calculates a refinement indicator
 * for them and creates the subspaces with high error indicators. Refinable subspaces that are not created
 * are stored for later reuse.
 */
class SubspaceGSGRefinement: public SubspaceRefinement {
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
	SubspaceGSGRefinement(AbstractRefinement* refinement, size_t dim);

	/**
	 * Refines a grid by adding additional subspaces according to a RefinementFunctor provided.
	 * Refines up to RefinementFunctor::getRefinementsNum() grid points if
	 * possible, and if their refinement value is larger than RefinementFunctor::start()
	 * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 */
	void freeRefineSubspace(GridStorage* storage,RefinementFunctor* functor);

protected:

	// availableGridPoints stores all refinable subspaces
	ErrorStorage availableSubspaces;
	// addedInLastRefinement vector to store the subspaces created in this refinement step.
	// based on these, we know where to look for for admissible subspaces in the next iterative refinement step
	ErrorVector addedInLastRefinement;
	// specifies how to look for refinable points
	bool firstRefinement;

	/**
	 * Checks a subspace for admissibility.
	 * A subspaces is called admissible, if all parent subspaces in all dimensions exist,
	 * i.e. if all parent in all dimensions exist for all grid points in a subspace.
	 *
	 * @param storage hashmap that stores the grid points
	 * @param subspace subspace to checked for admissibility.
	 */
	bool checkAdmissible(GridStorage* storage, index_type& subspace);


	/**
	 * Creates all subspaces that are passed in the maxErrorSubspaces array according to
	 * the specifications in the refinement functor.
	 *
	 * @param storage hashmap that stores the grid points
	 * @param errorStorage holds all refinable subspaces together with their error indicators.
	 * @param refinements_num number of subspaces to refine
	 * @param maxErrorSubspaces the array with the indices containing the level vectors
	 *  on which the subspaces should be created
	 * @param maxErrorValues the array with the corresponding indicator values
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 */
	virtual void refineSubspaceCollection(GridStorage* storage,
			ErrorStorage* errorStorage,
			ErrorVector* addedInLastStep,
			RefinementFunctor* functor);


	/**
	 * Examines the child subspaces of the subspaces from addedInLastRefinement and adds
	 * them to the error storage, if they are admissible (all grid points in the subspace have all parents in all dimensions).
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 * @param addedInLastRefinement a vector to store the subspaces created in this refinement step.
	 * @param errorStorage storage container for refinable subspaces, sorted by indicator value
	 */
	virtual void updateAdmissibleSubspaces(GridStorage* storage,
			RefinementFunctor* functor,
			ErrorVector* addedInLastRefinement,
			ErrorStorage* admissibleSubspaces);

};
} /* namespace base */
} /* namespace sg */
#endif /* SUBSPACEGSGREFINEMENT_HPP_ */
