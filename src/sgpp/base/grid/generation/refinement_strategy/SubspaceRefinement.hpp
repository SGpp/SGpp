/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef SUBSPACEREFINEMENT_HPP_
#define SUBSPACEREFINEMENT_HPP_


#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "dataStructures/ErrorContainer.hpp"
#include "dataStructures/ErrorStorage.hpp"
#include "RefinementDecorator.hpp"

namespace sg {
namespace base {


/**
 * Simple dimension adaptive refinement.
 * Selects hierarchical subspaces with at least one parent subspace, calculates a refinement indicator
 * for them and creates the subspaces with high error indicators.
 */
class SubspaceRefinement: public sg::base::RefinementDecorator {
public:

	typedef std::vector<index_type> SubspaceVector;
	typedef std::vector<index_type> GridPointVector;
	typedef std::vector<ErrorType> ErrorVector;


    /**
     * Constructor
     *
     * @param refinement object implementing the core functionality (e.g.
     * refinement with or without boundaries).
     */
	SubspaceRefinement(AbstractRefinement* refinement):RefinementDecorator(refinement){};

	/**
	 * Refines a grid by adding additional Subspaces according to a RefinementFunctor provided.
	 * Refines up to RefinementFunctor::getRefinementsNum() grid points if
	 * possible, and if their refinement value is larger than RefinementFunctor::start()
	 * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 */
	virtual void freeRefineSubspace(GridStorage* storage, RefinementFunctor* functor);

	/**
	 * This method creates a new subspace in the grid. It checks if some parents or
	 * children are needed in other dimensions.
	 *
	 * @param storage hashmap that stores the gridpoints
	 * @param index - the index containing the level vector for the new subspace.
	 * the index vector of the object has to be set to 1 for all dimensions!.
	 */
	virtual void createSubspace(GridStorage* storage, index_type& index, bool isLeaf);


	/*
	 *	Method to create grid points of a hierarchical subspace and store them in a vector starting at a given point.
	 *	If all grid points should be created, the index vector of the starting point has to be set to 1 in all dimensions.
	 *	Calls the recursive createAllGridPointsOfSubspaceHelper method with the input argumensts on level 0.
	 *
	 *	@param gridPoints vector, where the newly created grid points will be stored in
	 *	@param gridObject the grid point from which to start creating all grid points -
	 *	the index vector has to be set to 1 in all dimensions for this method
	 *	to create all grid points of a subspace.
	 *
	 */
	void createAllGridPointsOfSubspace(index_type& subspace, GridPointVector* gridPoints);

protected:

	/**
	 * Examines the grid subspaces, finds which ones are refinable and adds
	 * the error indicators of all points which belong to the same subspace.
	 * It returns an unsroted array with the refinements_num subspaces with the highest error indicator. (not sorted)
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 * @param errorStorage hash map that stores all subspaces that can be refined.
	 */
	virtual void collectRefinableSubspaces(GridStorage* storage,
			RefinementFunctor* functor,
			HashErrorStorage* errorStorage);


	/*
	 * Examines subspaces from a HashErrorStorage, selects the subspaces with the highest error indicator
	 * and stores the selected subspaces into an ErrorVector.
	 *
	 * @param errorStorage hashmap that stores the refinable subspaces with their error indicators
	 * @param refinements_num the amount of highest error subspaces to select for creation
	 * @param maxSubspaces a vector where the subspaces with the highest error indicators will be stored in
	 */
	virtual void selectHighestErrorSubspaces(HashErrorStorage* errorStorage,
											 size_t refinements_num,
											 ErrorVector* maxSubspaces);
	/**
	 * Creates all subspaces that are passed in the maxErrorSubspaces array according to
	 * the specifications in the refinement functor.
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 * @param refinements_num number of points to refine
	 * @param maxErrorSubspaces the array with the indices containing the level vectors
	 *  on which the subspaces should be created
	 * @param maxErrorValues the array with the corresponding indicator values
	 */
	virtual void refineSubspaceCollection(GridStorage* storage,
			RefinementFunctor* functor,
			size_t refinements_num,
			ErrorVector* maxErrorSubspaces);

private:

	/*
	 *	Recursive subroutine to create grid points of a hierarchical subspace and store them in a vector starting at a given point.
	 *	If all grid points should be created, the index vector of the starting point has to be set to 1 in all dimensions.
	 *
	 *	@param gridPoints vector, where the newly created grid points will be stored in
	 *	@param gridObject the grid point from which to start creating all grid points -
	 *	the index vector has to be set to 1 in all dimensions for this method
	 *	to create all grid points of a subspace.
	 *	@param dim - starting dimension
	 *
	 */
	void createAllGridPointsOfSubspaceHelper(GridPointVector* gridPoints,index_type& gridObject,size_t dim);

};


/*
 * Helper struct for priority queue. Sorts error objects smallest error first.
 */
struct smallestErrorFirst
{

	/*
	 * compare first argument with second,
	 * returning true if the first argument is bigger then the second
	 *
	 * @param lhs first error object
	 * @param rhs second error object
	 * @return true if first error object is bigger then second error object
	 */
	bool operator() (const ErrorType& lhs, const ErrorType& rhs) const
	{
		return lhs>rhs;
	}
};

} /* namespace base */
} /* namespace sg */
#endif /* SUBSPACEREFINEMENT_HPP_ */
