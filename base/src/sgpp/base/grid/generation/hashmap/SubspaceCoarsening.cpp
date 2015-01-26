/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "SubspaceCoarsening.hpp"
#include <vector>
#include <list>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {


	void SubspaceCoarsening::free_coarsen(GridStorage* storage,CoarseningFunctor* functor, DataVector* alpha)
	{
		SubspaceCoarseningErrorStorage errorStorage;
		SubspaceCoarseningErrorContainer* minContribSubspaces = new SubspaceCoarseningErrorContainer[functor->getRemovementsNum()];

		collectSubspacesToCoarsen(storage,functor,&errorStorage);

		pickLowestContributionSubspaces(&errorStorage,functor->getRemovementsNum(),minContribSubspaces);


		//DEBUG:
//		std::cout << "picked ";
//		for(size_t i=0; i<functor->getRemovementsNum(); ++i)
//		{
//			index_type* index = storage->get(minContribSubspaces[i].gridPoints[0]);
//			std::cout << index->toString() << ", is coarsenable " <<
//					minContribSubspaces[i].isCoarsenable << " has " << minContribSubspaces[i].gridPoints.size() << " points, " <<
//					"contributing " << minContribSubspaces[i].error << std::endl;
//		}


		//delete gridPoints
		std::list<size_t> pointsToDelete;

		// iterate over minContribSubspaces. add all points to be deleted to a list.
		for (size_t i = 0; i < functor->getRemovementsNum(); ++i) {

			//move all points of subspace to list
			for(std::vector<size_t>::iterator pointsIter = minContribSubspaces[i].gridPoints.begin();
					pointsIter != minContribSubspaces[i].gridPoints.end(); pointsIter++
					)
			{
				pointsToDelete.push_back(*pointsIter);
			}
		}
		//delete gridPoints
		std::vector<size_t> remainingIndex = storage->deletePoints(pointsToDelete);

		//delete all surpluses belonging to points deleted by coarsening from alpha vector;
		alpha->restructure(remainingIndex);

		//clean up
		delete[] minContribSubspaces;
	}

	void SubspaceCoarsening::collectSubspacesToCoarsen (GridStorage* gridStorage, CoarseningFunctor* functor, SubspaceCoarseningErrorStorage* errorStorage)
	{

		SubspaceCoarseningErrorStorage::iterator errorIter = errorStorage->begin();

		for (GridStorage::grid_map_iterator gridIter = gridStorage->begin();
				gridIter != gridStorage->end();
				++gridIter)
		{
			//we do only coarsen inner points.
			if (!gridIter->first->isInnerPoint()) {
				continue;
			}

			//try to find subspace of point in storage
			index_type* gridPoint = gridIter->first;
			errorIter = errorStorage->find(gridPoint);

			if (errorIter == errorStorage->end()) {
				// not found. add

				//create a new error container and store all coarsening information.
				SubspaceCoarseningErrorContainer errorContainer;
				errorContainer.isCoarsenable = gridIter->first->isLeaf();
				errorContainer.error = (*functor)(gridStorage, gridIter->second);
				errorContainer.gridPoints.push_back(gridIter->second);

				//add to error Container.
				errorStorage->insert(std::make_pair(gridIter->first,errorContainer));
			}else{

				// found - test if the grid point is a leaf (else it can not be coarsened)
				// or the subspace itself is already impossible to coarsen.
				if(errorIter->second.isCoarsenable && gridPoint->isLeaf()){

					// add to existing error.
					errorIter->second.error = (*functor)(gridStorage, gridIter->second);
					errorIter->second.gridPoints.push_back(gridIter->second);

				}else{

					//subspace can not be coarsened
					errorIter->second.isCoarsenable = false;
				}
			}
		}

		// DEBUG: print subspaces to coarsen
//		std::cout << "Subspaces:" << std::endl;
//		for (SubspaceCoarseningErrorStorage::iterator it = errorStorage->begin();
//				it != errorStorage->end(); ++it) {
//			std::cout << ((index_type) it->first).toString() << ", is coarsenable " <<
//					it->second.isCoarsenable << " has " << it->second.gridPoints.size() << " points, " <<
//					"contributing " << it->second.error << std::endl;
//		}
	}

	void SubspaceCoarsening::cleanUpErrorStorage(SubspaceCoarseningErrorStorage* errorStorage)
	{
		SubspaceCoarseningErrorStorage::iterator storageIter;

		for(storageIter = errorStorage->begin(); storageIter != errorStorage->end(); ++storageIter)
		{
			//@TODO:check if index is not moved by removing!
			// subspace can not be coarsened. clean it up.
			if(!storageIter->second.isCoarsenable)
			{
				errorStorage->erase(storageIter);
			}
		}

	}

	void SubspaceCoarsening::pickLowestContributionSubspaces(SubspaceCoarseningErrorStorage* errorStorage,
										 size_t removementsNum,
										 SubspaceCoarseningErrorContainer* minContribSubspaces )
	{
		//initialize
		size_t maxIndex = 0;
		//create and init helper array
		for (size_t i = 0; i < removementsNum; ++i) {
			minContribSubspaces[i] = SubspaceCoarseningErrorContainer();
			minContribSubspaces[i].error = INFINITY;
		}

		//iterate over all subspaces
		SubspaceCoarseningErrorStorage::iterator storageIter;
		for(storageIter = errorStorage->begin(); storageIter != errorStorage->end(); ++storageIter)
		{
			//if subspace can be coarsened
			//and if its contribution is smaller then one of the others, insert it
			if (storageIter->second.isCoarsenable &&
					storageIter->second.error < minContribSubspaces[maxIndex].error) {

				//override the entry with the highest contribution to the grid
				minContribSubspaces[maxIndex] = storageIter->second;

				//find a new maximum
				maxIndex = getMaxErrorElem(minContribSubspaces,removementsNum);
			}
		}
	}

	size_t SubspaceCoarsening::getMaxErrorElem(SubspaceCoarseningErrorContainer* minContribSubspaces, size_t removementsNum)
	{
		size_t maxIndex = 0;

		for (size_t i = 1; i < removementsNum; i++) {
			if (minContribSubspaces[i].error > minContribSubspaces[maxIndex].error) {
				maxIndex = i;
			}
		}

		return maxIndex;
	}


} /* namespace base */
} /* namespace SGPP */
