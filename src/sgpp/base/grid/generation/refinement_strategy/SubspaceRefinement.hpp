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


class SubspaceRefinement: public sg::base::RefinementDecorator {
public:

	typedef std::vector<index_type> SubspaceVector;
	typedef std::vector<index_type> GridPointVector;
	typedef std::vector<ErrorType> ErrorVector;


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

	void createAllGridPointsOfSubspace(index_type& subspace, GridPointVector* gridPoints);

protected:

	/**
	 * Examines the grid points, finds which ones are refinable and adds
	 * the error indicators of all points which belong to the same subspace.
	 * It returns an unsroted array with the refinements_num subspaces with the highest error indicator. (not sorted)
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 * @param refinements_num number of points to refine
	 * @param maxErrorSubspaces the array of hash grid indices, containing the level vector of the subspace
	 * @param maxErrorValues the array with the indicator values corresponding to the level of the subspace.
	 */
	virtual void collectRefinableSubspaces(GridStorage* storage,
			RefinementFunctor* functor,
			HashErrorStorage* errorStorage);


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

	/*
	 * calculates the
	 */
	//bool getParentLevelAndIndex(index_type* child,size_t dim);


	/**
	 * Sets all the elements of the index vector of a grid point to 1.
	 * @param gridPoint pointer to the grid point with the index array to be changed
	 */
	//void resetIndexVector(index_type* gridPoint);

private:

	void createAllGridPointsOfSubspaceHelper(GridPointVector* gridPoints,index_type& gridObject,size_t dim);

	/**
	 * recursive function to create all points on a subspace.
	 *
	 * @param storage
	 * @param index
	 * @param dim
	 *
	 */

	//void createGridPointWithParents(GridStorage* storage, index_type& child);

};

struct smallestErrorFirst
{
	bool operator() (const ErrorType& lhs, const ErrorType& rhs) const
	{
		return lhs>rhs;
	}
};

} /* namespace base */
} /* namespace sg */
#endif /* SUBSPACEREFINEMENT_HPP_ */
