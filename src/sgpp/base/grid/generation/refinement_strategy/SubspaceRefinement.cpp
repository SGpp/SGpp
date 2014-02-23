/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "SubspaceRefinement.hpp"
#include <queue>

namespace sg {
namespace base {

void SubspaceRefinement::collectRefinableSubspaces(GridStorage* storage,
											   RefinementFunctor* functor,
											   HashErrorStorage* subspaceError)
{
	//work through all refinable points
	index_type index;
	GridStorage::grid_map_iterator end_iter = storage->end();

	// start iterating over whole grid
	for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
		index = *(iter->first);

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


			if (child_iter == end_iter) {
				//grid point is not contained in the grid => its subspace is not part of the grid
				RefinementFunctor::value_type gridPointError = (*functor)(storage, iter->second);

				//check if the grid point's subspace is already in map
				ErrorType errorContainer(gridPointError,index);
				//therefore first copy and reset index vector to search for subspace
				errorContainer.resetIndexVector();
				HashErrorStorage::grid_map_iterator subSpaceErrorIterator = subspaceError->find(&errorContainer);

				//not found
				if(subSpaceErrorIterator == subspaceError->end())
				{
					//insert new subspace
					subspaceError->insert(errorContainer);
				}
				//found
				else
				{
					//add error to subspace
					*(subSpaceErrorIterator->first) += gridPointError;
				}
			}

			// test existance of right child
			index.set(d, source_level + 1, 2 * source_index + 1);
			child_iter = storage->find(&index);

			if (child_iter == end_iter) {
				RefinementFunctor::value_type gridPointError = (*functor)(storage, iter->second);

				//check if the grid point's subspace is already in map
				ErrorType errorContainer(gridPointError,index);
				//therefore first copy and reset index vector to search for subspace
				errorContainer.resetIndexVector();
				HashErrorStorage::grid_map_iterator subSpaceErrorIterator = subspaceError->find(&errorContainer);

				//not found
				if(subSpaceErrorIterator == subspaceError->end())
				{
					//insert new subspace
					subspaceError->insert(errorContainer);
				}
				//found
				else
				{
					//add error to subspace
					(*subSpaceErrorIterator->first) += gridPointError;
				}
			}
			// reset current grid point in dimension d
			index.set(d, source_level, source_index);
		}
	}



//	//DEBUG print map
//	std::cout << "Map contains:\n";
//	for (HashErrorStorage::grid_map_iterator errorIter = subspaceError->begin();errorIter != subspaceError->end(); errorIter++)
//	{
//		std::cout << "Subspace" << ((ErrorType) (errorIter->first)).toString() << std::endl;
//
//	}
}

void SubspaceRefinement::refineSubspaceCollection(GridStorage* storage, RefinementFunctor* functor, size_t refinements_num, ErrorVector* maxErrorSubspaces) {

	RefinementFunctor::value_type maxErrorValue;
	index_type maxErrorSubspaceIndex;

	double threshold = functor->getRefinementThreshold();

	// now refine all grid points which satisfy the refinement criteria
	for (ErrorVector::iterator iter = maxErrorSubspaces->begin(); iter != maxErrorSubspaces->end(); ++iter) {
		maxErrorValue = (*iter).getContribPerPoint();
		maxErrorSubspaceIndex = (*iter);

		if (maxErrorValue > functor->start() && fabs(maxErrorValue) >= threshold) {
			createSubspace(storage,maxErrorSubspaceIndex,false);
		}
	}
}

void SubspaceRefinement::freeRefineSubspace(GridStorage* storage, RefinementFunctor* functor) {

	//nothing there => nothing to refine
	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}

	// the functor->getRefinementsNum() largest subspaces should be refined.
	size_t refinements_num = functor->getRefinementsNum();
	//store the subspaces here with the largest error indicator.
	ErrorVector maxSubspaces;

	//accumulate error on refinable subspaces
	HashErrorStorage errorStorage(storage->dim());
	collectRefinableSubspaces(storage,functor,&errorStorage);

	//select refinements_num highest indicator subspaces
	selectHighestErrorSubspaces(&errorStorage,refinements_num,&maxSubspaces);

	//refine all subspaces which satisfy the refinement criteria
	refineSubspaceCollection(storage,functor,refinements_num,&maxSubspaces);

}

void SubspaceRefinement::createSubspace(GridStorage* storage, index_type& gridIndex, bool isLeaf) {



	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}


	index_type parentSubspace;
	GridStorage::grid_map_iterator storageIndex;
	//first see if all parent subspaces are available
	for(size_t dim = 0; dim < gridIndex.dim(); ++dim){

		level_t parentLevel;
		index_t parentIndex;
		gridIndex.getParentLevelAndIndex(&parentLevel,&parentIndex,dim);
		parentSubspace = gridIndex;
		parentSubspace.set(dim,parentLevel,parentIndex);

		storageIndex = storage->find(&parentSubspace);

		if ( parentLevel>0 && storageIndex == storage->end())
		{
			//parent in dimension dim not found - create the parent subspace
			createSubspace(storage,parentSubspace,true);

		}else if (parentLevel>0 && storageIndex->first->isLeaf()) {
			//remove the leaf property from all future parents

			GridPointVector removeLeafProperty;
			index_type newSubspace = storageIndex->first;
			newSubspace.resetIndexVector();
			createAllGridPointsOfSubspace(newSubspace,&removeLeafProperty);

			for(GridPointVector::iterator gridPointIter = removeLeafProperty.begin(); gridPointIter != removeLeafProperty.end(); ++gridPointIter)
			{
				storageIndex = storage->find(&(*gridPointIter));
				storageIndex->first->setLeaf(false);
			}
		}
	}

	//get all indices we want to create.
	GridPointVector gridPointsToCreate;
	createAllGridPointsOfSubspace(gridIndex,&gridPointsToCreate);

	//create all the grid points
	for(GridPointVector::iterator gridPointIter = gridPointsToCreate.begin(); gridPointIter != gridPointsToCreate.end(); ++gridPointIter)
	{
		createGridpoint(storage,*gridPointIter);
		GridStorage::grid_map_iterator storageIndex = (storage->find(&(*gridPointIter)));
		storageIndex->first->setLeaf(isLeaf);
	}

}


void SubspaceRefinement::selectHighestErrorSubspaces(HashErrorStorage* subspaceError,
		size_t refinements_num,
		ErrorVector* maxErrorSubspaces)
{
	//sort by adding indices to a limited size priority queue, sorted  smallest error subspaces first.
	std::priority_queue<ErrorType,ErrorVector,smallestErrorFirst> largestErrorObjects;
	size_t size = 0;


	//iterate over all refinable subspaces
	for (HashErrorStorage::grid_map_iterator subspaceErrorIterator = subspaceError->begin(); subspaceErrorIterator != subspaceError->end(); subspaceErrorIterator++)
	{
		// simply add subspace to priority queue, if the priority queue has less elements then refinable points.
		if(size < refinements_num)
		{
		  largestErrorObjects.push(*(subspaceErrorIterator->first));
		}

		//if the priority queue is full and the error indicator of the current subspace is larger
		//then the the error indicator for the subspace with smallest error indicator in the queue,
		//replace it.
		if( size >= refinements_num && (largestErrorObjects.top() < *(subspaceErrorIterator->first)) )
		{
			ErrorType error = largestErrorObjects.top();
			largestErrorObjects.pop();
			largestErrorObjects.push(*(subspaceErrorIterator->first));
		}

		++size;
	}

	//transfer indices to a vector
	while(!largestErrorObjects.empty())
	{
		maxErrorSubspaces->push_back(largestErrorObjects.top());
		largestErrorObjects.pop();
	}

//	//DEBUG print errors
//	std::cout << "Selected:\n";
//	for (ErrorVector::iterator iter = maxErrorSubspaces->begin(); iter != maxErrorSubspaces->end();++iter) {
//		ErrorType error = *iter;
//		std::cout << "Subspace" << error.toString() << std::endl;
//	}
}

void SubspaceRefinement::createAllGridPointsOfSubspace(index_type& subspace,
		GridPointVector* gridPoints) {

	createAllGridPointsOfSubspaceHelper(gridPoints,subspace,0);
}

void SubspaceRefinement::createAllGridPointsOfSubspaceHelper(
		GridPointVector* gridPoints, index_type& storageIndex, size_t dim)
{

	//walk through all dimensions
	if (dim < storageIndex.dim()) {

		//get level of subspace on that dimension
		index_t index;
		level_t level;
		storageIndex.get(dim,level,index);
		index = 1;

		//iterate over all allowed indices on that level, in that dimension
		while(index < static_cast <size_t>( 1 << level)){

			//set gridpoint's index accordingly
			storageIndex.set(dim,level,index);
			//recursive call, so that we iterate over all indices on all levels in all dimensions
			createAllGridPointsOfSubspaceHelper(gridPoints,storageIndex,dim+1);
			//move to next admissible
			index=index+2;
		}

	} else {
		//reached end of recursion scheme.
		//Add grid point to grid points contained in the subspace.
		gridPoints->push_back(storageIndex);
	}
}

} /* namespace base */
} /* namespace sg */
