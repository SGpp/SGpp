/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "PredictiveSubspaceGSGRefinement.hpp"

using namespace std;

namespace sg {
namespace base {

PredictiveSubspaceGSGRefinement::PredictiveSubspaceGSGRefinement(RefinementDecorator* decorator):SubspaceGSGRefinement(decorator)
{}

void sg::base::PredictiveSubspaceGSGRefinement::freeRefineSubspace(
		GridStorage* storage, PredictiveRefinementIndicator* errorIndicator)
{
	//nothing there => nothing to refine
	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}

	// the functor->getRefinementsNum() largest subspaces should be refined.
	size_t refinements_num = errorIndicator->getRefinementsNum();
	//subspaces
	index_type* maxSubspaces = new index_type[refinements_num];
	//subspace errors
	RefinementFunctor::value_type* errorsPerSubspace = new RefinementFunctor::value_type[refinements_num];

	//accumulate error on refinable subspaces
	if (availableSubspaces.empty()) {
		collectRefinableSubspaces(storage,errorIndicator,&availableSubspaces);

		//set all refinable subspaces as admissible, where parents in all dimensions exist
		for(SubspaceErrorStorage::iterator errorIter = availableSubspaces.begin();
				errorIter != availableSubspaces.end(); ++errorIter)
		{
			index_type subspace = errorIter->first;
			errorIter->second.setAdmissible(checkAdmissible(storage,subspace));
		}

	}else{

		// check if refinement in last step made nonadmissible Subspaces admissible
		for(SubspaceVector::iterator lastAddedIter = addedInLastRefinement.begin();
				lastAddedIter !=  addedInLastRefinement.end();
				++lastAddedIter)
		{
			index_type neighbour = *lastAddedIter;
			index_t index = 1;
			level_t level = 1;

			for (size_t dim = 0; dim < storage->dim(); ++dim)
			{
				neighbour.get(dim,level,index);
				neighbour.set(dim,level+1,index);
				SubspaceErrorStorage::iterator errorIter = availableSubspaces.find(neighbour);
				errorIter->second.setAdmissible(checkAdmissible(storage,neighbour));
				neighbour = *lastAddedIter;
			}
		}

		updateAdmissibleSubspaces(storage,errorIndicator,&addedInLastRefinement,&availableSubspaces);
	}

	//DEBUG: print all elements

	std::cout << "\n\navailable subspaces" << "\n";
	std::cout << "=================================================" << "\n";
	for(SubspaceErrorStorage::iterator errorIter = availableSubspaces.begin(); errorIter != availableSubspaces.end(); ++errorIter)
	{
		std::cout << ((index_type) errorIter->first).toString() << " , " << errorIter->second.toString() << "\n";
	}

	//select refinements_num highest indicator subspaces
	selectHighestErrorSubspaces(&availableSubspaces,refinements_num,maxSubspaces,errorsPerSubspace);

	//DEBUG: selected
	std::cout << "\nselected subspaces" << "\n";
	std::cout << "=================================================" << "\n";
	for (size_t i = 0; i < refinements_num; ++i)
	{
		std::cout << ((index_type) maxSubspaces[i]).toString() << " , " << errorsPerSubspace[i] << "\n";
	}


	//refine all subspaces which satisfy the refinement criteria
	refineSubspaceCollection(storage,&availableSubspaces,&addedInLastRefinement,errorIndicator,refinements_num,maxSubspaces,errorsPerSubspace);

	delete [] maxSubspaces;
	delete [] errorsPerSubspace;


}

void sg::base::PredictiveSubspaceGSGRefinement::updateAdmissibleSubspaces(
		GridStorage* storage, PredictiveRefinementIndicator* errorIndicator,
		SubspaceVector* addedInLastRefinement,
		SubspaceErrorStorage* availableSubspaces)
{
	//	std::cout << "updating admissible subspaces \n";

		SubspaceVector newSubspaces;
		//go through all subspaces refined in last refinement step.
		for (SubspaceVector::iterator addedInLastIter = addedInLastRefinement->begin();
				addedInLastIter != addedInLastRefinement->end(); ++addedInLastIter)
		{

			//create all gridPoints in this subspace.
			index_type subspace = *addedInLastIter;
			GridPointVector gridPoints;
			createAllGridPointsOfSubspace(subspace,&gridPoints);

			//for each grid point, go through all dimensions
			for(GridPointVector::iterator pointIter = gridPoints.begin(); pointIter != gridPoints.end();++pointIter)
			{
				//find out if children exist in
				// 1) a refinable, non admissible subspace,
				// 2) a non existing subspace

				index_t index = 1;
				level_t level = 1;

				for(size_t dim = 0; dim < storage->dim(); ++dim)
				{
					//find neighbouring subspace.
					index_type neighbourSubspace = *pointIter;
					(*pointIter).get(dim,level,index);

					neighbourSubspace.set(dim,level+1,1);

					SubspaceErrorStorage::iterator neighbourIter = availableSubspaces->find(neighbourSubspace);

					if(neighbourIter == availableSubspaces->end())
					{

	//					std::cout << "did not find " << neighbourSubspace.toString() << "\n";
						//if the neighbouring subspace does not exist, we will create it.
						resetIndexVector(&neighbourSubspace);
						std::pair<SubspaceErrorStorage::iterator,bool> insertionLocation =
								availableSubspaces->insert(std::make_pair(neighbourSubspace,ErrorContainer()));
						newSubspaces.push_back(neighbourSubspace);
						neighbourIter = insertionLocation.first;
					}


					// add error left child;

					// get levels and indices on level and dim;
					(*pointIter).get(dim,level,index);
					index_type child = (*pointIter);
					child.set(dim, level + 1, 2 * index - 1);

					(neighbourIter->second)+= (*errorIndicator)(&child);

					// add error right child

					// test existance of right child
					(*pointIter).get(dim,level,index);
					child = (*pointIter);
					child.set(dim, level + 1, 2 * index + 1);

					(neighbourIter->second)+= (*errorIndicator)(&child);
					//std::cout << "updated error to " <<  neighbourIter->second.toString() << "\n";
				}
			}
		}

		//check added subspaces for admissibility.
		std::cout << "checking subspaces for admissibility \n";
		for (SubspaceVector::iterator newSubspaceIter = newSubspaces.begin();
					newSubspaceIter != newSubspaces.end(); ++newSubspaceIter)
		{
			SubspaceErrorStorage::iterator errorIter = availableSubspaces->find(*newSubspaceIter);
			index_type subspace = errorIter->first;
			errorIter->second.setAdmissible(checkAdmissible(storage, subspace));
			std::cout << ((index_type) errorIter->first).toString() << " ; " <<  errorIter->second.toString() << "\n";
		}


		// DEBUG:: print all new subspaces

}

void sg::base::PredictiveSubspaceGSGRefinement::collectRefinableSubspaces(
		GridStorage* storage, PredictiveRefinementIndicator* errorIndicator,
		SubspaceErrorStorage* subspaceError)
{
	//storage for accumulated errors on each subspace
	//SubspaceErrorStorage subspaceError;


	//work through all refinable points

	index_type index;
	GridStorage::grid_map_iterator end_iter = storage->end();

	// start iterating over whole grid
	for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
		index = *(iter->first);

		cout<< "analyzing point"<< index.toString() << "\n";
		GridStorage::grid_map_iterator child_iter;

		// @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
		// check for each grid point whether it can be refined (i.e., whether not all kids exist yet)
		// if yes, check whether it belongs to the refinements_num largest ones
		for (size_t d = 0; d < storage->dim(); d++) {
			index_t source_index;
			level_t source_level;
			index.get(d, source_level, source_index);

			// test existence of left child
			index.set(d, source_level + 1, 2 * source_index - 1);
			child_iter = storage->find(&index);

			//cout <<"testing existence of left child " << index.toString() <<"\n";

			// if there no more grid points --> test if we should refine the grid
			if (child_iter == end_iter) {
				RefinementFunctor::value_type gridPointError = (*errorIndicator)(&index);

				//check if already in map
				//cout << "running find operation \n";
				SubspaceErrorStorage::iterator subSpaceErrorIterator = subspaceError->find(index);

				//not found
				if(subSpaceErrorIterator == subspaceError->end())
				{
					//cout << "inserting new subspace " << index.toString() <<"\n";
					//insert new subspace
					index_type indexCopy = index;
					//reset the index vector to one.
					resetIndexVector(&indexCopy);
					//cout << indexCopy.toString() << "\n";
					subspaceError->insert(make_pair(indexCopy, ErrorContainer(gridPointError)));
				}
				//found
				else
				{

					//cout << "adding to old subspace" <<  ((index_type)(subSpaceErrorIterator->first)).toString() << "\n" ;
					//add error to subspace
					//cout << "old " << ((errorContainer) subSpaceErrorIterator->second).toString();
					subSpaceErrorIterator->second += gridPointError;
					//cout << " new " << ((errorContainer) subSpaceErrorIterator->second).toString() << "\n";
				}
				//break;
			}

			// test existance of right child
			index.set(d, source_level + 1, 2 * source_index + 1);
			//cout <<"testing existence of right child " << index.toString() <<"\n";
			child_iter = storage->find(&index);

			if (child_iter == end_iter) {
				RefinementFunctor::value_type gridPointError = (*errorIndicator)(&index);

				//check if already in map
				SubspaceErrorStorage::iterator subSpaceErrorIterator = subspaceError->find(index);

				//not found
				if(subSpaceErrorIterator == subspaceError->end())
				{
					//cout << "inserting new subspace " << index.toString() <<"\n";
					//insert new subspace
					index_type indexCopy = index;
					//reset the index vector to one.
					resetIndexVector(&indexCopy);
					//cout << indexCopy.toString() << "\n";
					subspaceError->insert(make_pair(indexCopy, ErrorContainer(gridPointError)));
				}
				//found
				else
				{
					//cout << "adding to old subspace" <<  ((index_type)(subSpaceErrorIterator->first)).toString() << "\n" ;
					//add error to subspace
					//cout << "old " << ((errorContainer) subSpaceErrorIterator->second).toString();
					subSpaceErrorIterator->second += gridPointError;
					//cout << " new " << ((errorContainer) subSpaceErrorIterator->second).toString() << "\n";
				}
				//break;
			}
			// reset current grid point in dimension d
			index.set(d, source_level, source_index);
		}
	}

	cout << "Map contains:\n";

	//debug method - print map
	for (SubspaceErrorStorage::iterator errorIter = subspaceError->begin();errorIter != subspaceError->end(); errorIter++)
	{
		cout << "Subspace" << ((index_type) (errorIter->first)).toString() << " with error " <<  ((ErrorContainer) errorIter->second).getContribPerPoint() <<"\n";

	}
}


void sg::base::PredictiveSubspaceGSGRefinement::refineSubspaceCollection(GridStorage* storage,
		SubspaceErrorStorage* errorStorage,
		SubspaceVector* addedInLastStep,
		PredictiveRefinementIndicator* errorIndicator,
		size_t refinements_num,
		index_type* maxErrorSubspaces,
		RefinementFunctor::value_type* maxErrorValues)
{
	//cout << "refining subspace collection\n";
	RefinementFunctor::value_type maxErrorValue;
	index_type maxErrorSubspaceIndex;
	// now refine all grid points which satisfy the refinement criteria
	double threshold = errorIndicator->getRefinementThreshold();

	for (size_t i = 0; i < refinements_num; i++) {
		maxErrorValue = maxErrorValues[i];
		maxErrorSubspaceIndex = maxErrorSubspaces[i];

		if (maxErrorValue > errorIndicator->start() && fabs(maxErrorValue) >= threshold) {
			createSubspace(storage,maxErrorSubspaceIndex);

			//insert the subspaces that satisfy the refinement criteria into the storage
			// and remove newly created subspaces from available subspaces.
			errorStorage->erase(maxErrorSubspaceIndex);
			addedInLastStep->push_back(maxErrorSubspaceIndex);
		}
	}
}


} /* namespace base */
} /* namespace sg */
