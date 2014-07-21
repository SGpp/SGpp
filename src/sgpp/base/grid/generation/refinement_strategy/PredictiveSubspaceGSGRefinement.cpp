/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "PredictiveSubspaceGSGRefinement.hpp"

#include <typeinfo>

namespace sg {
namespace base {

void PredictiveSubspaceGSGRefinement::updateAdmissibleSubspaces(GridStorage* storage,
		RefinementFunctor* functor,
		ErrorVector* addedInLastRefinement,
		ErrorStorage* availableSubspaces)
{

	//perform dynamic down cast to see, if the functor is a valid predictive refinement indicator;

	PredictiveRefinementIndicator* errorIndicator = dynamic_cast<PredictiveRefinementIndicator*>(functor);

	HashErrorStorage* errorStorage = availableSubspaces->getHashErrorStorage();
	HashErrorStorage newSubspaces(storage->dim());

	//go through all subspaces refined in last refinement step.
	for (ErrorVector::iterator addedInLastIter = addedInLastRefinement->begin();
			addedInLastIter != addedInLastRefinement->end(); ++addedInLastIter)
	{
		//for each of them, find the neighbouring subspaces
		index_t index = 1;
		level_t level = 1;

		for(size_t dim = 0; dim < storage->dim(); ++dim)
		{
			//set to neighbouring subspace's level in dimension dim
			ErrorType neighbourSubspace = *addedInLastIter;
			neighbourSubspace.get(dim,level,index);
			neighbourSubspace.set(dim,level+1,index);

			//find out if the neighbour is already available.
			HashErrorStorage::grid_map_iterator newSubspacesIter = newSubspaces.find(&neighbourSubspace);
			HashErrorStorage::grid_map_iterator errorStorageIter = errorStorage->find(&neighbourSubspace);
			if(newSubspacesIter == newSubspaces.end() && errorStorageIter == errorStorage->end() && checkAdmissibility(storage,neighbourSubspace))
			{
				//not found. insert that subspace.
				GridPointVector gridPoints;
				ErrorType helper = neighbourSubspace;
				createAllGridPointsOfSubspace(helper,&gridPoints);
				//here the error contributions will be stored.
				ErrorType errorContribution;


				for(GridPointVector::iterator pointIter = gridPoints.begin(); pointIter != gridPoints.end();++pointIter)
				{
					// functor in each dimension is equal -> error Contribution of each subspace to each neighbouring subspace is equal.
					RefinementFunctor::value_type error =  (*errorIndicator)(&(*pointIter));

					//add error contrib for left child and for right child
					//careful! += 2*error is wrong, it does not increase the internal contributions counter!
					errorContribution += error;
				}

				neighbourSubspace+=errorContribution;
				newSubspaces.insert(neighbourSubspace);
			}
		}

		//update available subspaces
		availableSubspaces->updateErrors(&newSubspaces);
	}
}

void PredictiveSubspaceGSGRefinement::collectRefinableSubspaces(GridStorage* storage,
											   RefinementFunctor* functor,
											   HashErrorStorage* subspaceError)
{

	//perform dynamic down cast to see, if the functor is a valid predictive refinement indicator;
	PredictiveRefinementIndicator* errorIndicator = dynamic_cast<PredictiveRefinementIndicator*>(functor);

	//work through all refinable points

	index_type index;
	GridStorage::grid_map_iterator end_iter = storage->end();

	// start iterating over whole grid
	for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
		index = *(iter->first);

		GridStorage::grid_map_iterator neighbourIter;

		// @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
		// check for each grid point whether it can be refined (i.e., whether not all kids exist yet)
		// if yes, check whether it belongs to the refinements_num largest ones
		for (size_t d = 0; d < storage->dim(); d++) {
			index_t source_index;
			level_t source_level;
			index.get(d, source_level, source_index);

			// look for neighbour subspace
			index.set(d, source_level + 1,1);
			neighbourIter = storage->find(&index);


			// if there no more grid points --> test if we should refine the grid
			if (neighbourIter == end_iter) {
				//check if already in map
				//therefore first copy and reset index vector to search for subspace
				ErrorType errorContainer(index);
				errorContainer.resetIndexVector();
				HashErrorStorage::grid_map_iterator subSpaceErrorIterator = subspaceError->find(&errorContainer);

				if(subSpaceErrorIterator == subspaceError->end() && checkAdmissibility(storage,errorContainer))
				{
					//not found -> insert new subspace
					subspaceError->insert(errorContainer);
				}
			}

			// reset current grid point in dimension d
			index.set(d, source_level, source_index);
		}
	}

	//calculate error
	for (HashErrorStorage::grid_map_iterator errorIter = subspaceError->begin();errorIter != subspaceError->end(); ++errorIter)
	{

		GridPointVector gridPoints;
		ErrorType helper = errorIter->first;
		createAllGridPointsOfSubspace(helper,&gridPoints);
		//here the error contributions will be stored.
		ErrorType errorContribution = *(errorIter->first);

		for(GridPointVector::iterator pointIter = gridPoints.begin(); pointIter != gridPoints.end();++pointIter)
		{
			RefinementFunctor::value_type error =  (*errorIndicator)(&(*pointIter));
			errorContribution += error;
		}
		//update the error in this subspace. this invalidates the iterator.
		subspaceError->update(errorContribution,errorIter->second);
		errorIter = subspaceError->find(&errorContribution);
	}


//	//DEBUG method - print map
//	std::cout << "Map contains:\n";
//	for (HashErrorStorage::grid_map_iterator errorIter = subspaceError->begin();errorIter != subspaceError->end(); errorIter++)
//	{
//		std::cout << "Subspace" << ((ErrorType) (errorIter->first)).toString() << std::endl;
//	}
}

} /* namespace base */
} /* namespace sg */
