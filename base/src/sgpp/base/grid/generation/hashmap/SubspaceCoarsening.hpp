/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef SUBSPACECOARSENING_HPP_
#define SUBSPACECOARSENING_HPP_

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/dataStructures/SubspaceCoarseningErrorContainer.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>


namespace sg {
namespace base {

/**
 * Dimension adaptive coarsening class, which removes complete subspaces of a sparse grid, only
 * inner grid points can be removed
 */
class SubspaceCoarsening: public HashCoarsening {

public:

    /**
     * Performs coarsening on grid. It's possible to remove a certain number
     * of subspaces in one coarsening step. This number is specified within the
     * declaration of the coarsening functor. Also the coarsening threshold is
     * specified in the coarsening functor. ONLY INNER GRID POINTS WILL
     * BE REMOVED!
     *
     * @param storage hashmap that stores the grid points
     * @param functor a function used to determine if refinement is needed
     * @param alpha pointer to the gridpoints' coefficients removed points must also be considered in this vector
     */
	void free_coarsen(GridStorage* storage,CoarseningFunctor* functor, DataVector* alpha);

protected:


	/**
	 * Examines the subspaces, finds which ones which can be removed and adds
	 * the error indicators of all points which belong to the same subspace.
	 * It returns an unsroted array with the refinements_num subspaces with the smallest error indicator.
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 * @param errorStorage tree that stores all subspaces that can be removed.
	 */
	void collectSubspacesToCoarsen (GridStorage* gridStorage, CoarseningFunctor* functor, SubspaceCoarseningErrorStorage* errorStorage);

	/*
	 * Examines subspaces from a SubspaceCoarseningErrorStorage, selects the subspaces with the smallest error indicator
	 * and stores the selected subspaces into an array of SubspaceCoarseningErrorContainers.
	 *
	 * @param errorStorage hashmap that stores the refinable subspaces with their error indicators
	 * @param refinements_num the amount of highest error subspaces to select for creation
	 * @param maxSubspaces a vector where the subspaces with the highest error indicators will be stored in
	 */
	void pickLowestContributionSubspaces(SubspaceCoarseningErrorStorage* errorStorage,
			size_t removementsNum,
			SubspaceCoarseningErrorContainer* minContribSubspaces );

	/**
	 * Removes all subspaces that can not be coarsened from the error storage.
	 *
	 * @param errorStorage storage of the error
	 */
	void cleanUpErrorStorage(SubspaceCoarseningErrorStorage* errorStorage);

private:


	/*
	 *	searches an array of subspaces for the subspace with the highest error.
	 *
	 *	@param minContribSubspaces array of subspaces to select from
	 *	@param removementsNum size of the array
	 *	@return position of the subspace with the highest error indicator in the array
	 *
	 */
	size_t getMaxErrorElem(SubspaceCoarseningErrorContainer* minContribSubspaces, size_t removementsNum);
};

} /* namespace base */
} /* namespace sg */
#endif /* SUBSPACECOARSENING_HPP_ */
