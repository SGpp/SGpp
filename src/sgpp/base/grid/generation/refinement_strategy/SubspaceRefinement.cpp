/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "SubspaceRefinement.hpp"
#include <queue>

using namespace std;


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

//		std::cout<< "analyzing point"<< index.toString() << "\n";
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

//			std::cout <<"testing existence of left child " << index.toString() <<"\n";

			// if there no more grid points --> test if we should refine the grid
			if (child_iter == end_iter) {
				RefinementFunctor::value_type gridPointError = (*functor)(storage, iter->second);

				//check if already in map
				//therefore first copy and reset index vector to search for subspace
				ErrorType errorContainer(gridPointError,index);
				errorContainer.resetIndexVector();
				HashErrorStorage::grid_map_iterator subSpaceErrorIterator = subspaceError->find(&errorContainer);

				//not found
				if(subSpaceErrorIterator == subspaceError->end())
				{
//					std::cout << "inserting new subspace " << errorContainer.toString() <<"\n";
					//insert new subspace
//					std::cout << "here:" << errorContainer.toString() << "\n";
					subspaceError->insert(errorContainer);
				}
				//found
				else
				{
//					std::cout << "adding to old subspace" <<  ((ErrorType)(subSpaceErrorIterator->first)).toString() << "\n" ;
					//add error to subspace
//					std::cout << "old " << ((ErrorType) subSpaceErrorIterator->first).toString();
					*(subSpaceErrorIterator->first) += gridPointError;
//					std::cout << " new " << ((ErrorType) subSpaceErrorIterator->first).toString() << "\n";
				}
				//break;
			}

			// test existance of right child
			index.set(d, source_level + 1, 2 * source_index + 1);
//			cout <<"testing existence of right child " << index.toString() <<"\n";
			child_iter = storage->find(&index);

			if (child_iter == end_iter) {
				RefinementFunctor::value_type gridPointError = (*functor)(storage, iter->second);

				//check if already in map
				//therefore first copy and reset index vector to search for subspace
				ErrorType errorContainer(gridPointError,index);
				errorContainer.resetIndexVector();
				HashErrorStorage::grid_map_iterator subSpaceErrorIterator = subspaceError->find(&errorContainer);

				//not found
				if(subSpaceErrorIterator == subspaceError->end())
				{
//					std::cout << "inserting new subspace " << errorContainer.toString() <<"\n";
					//insert new subspace
					//cout << indexCopy.toString() << "\n";
					subspaceError->insert(errorContainer);
//					cout << errorContainer.toString() << "\n";
				}
				//found
				else
				{
//					cout << "adding to old subspace" <<  ((index_type)(subSpaceErrorIterator->first)).toString() << "\n" ;
					//add error to subspace
//					cout << "old " << ((ErrorType) subSpaceErrorIterator->first).toString();
					(*subSpaceErrorIterator->first) += gridPointError;
//					cout << " new " << ((ErrorType) subSpaceErrorIterator->first).toString() << "\n";
				}
				//break;
			}
			// reset current grid point in dimension d
			index.set(d, source_level, source_index);
		}
	}

	cout << "Map contains:\n";

	//debug method - print map
	for (HashErrorStorage::grid_map_iterator errorIter = subspaceError->begin();errorIter != subspaceError->end(); errorIter++)
	{
		cout << "Subspace" << ((ErrorType) (errorIter->first)).toString() << std::endl;

	}
}

void SubspaceRefinement::refineSubspaceCollection(GridStorage* storage, RefinementFunctor* functor, size_t refinements_num, ErrorVector* maxErrorSubspaces) {

	//cout << "refining subspace collection\n";
	RefinementFunctor::value_type maxErrorValue;
	index_type maxErrorSubspaceIndex;
	// now refine all grid points which satisfy the refinement criteria
	double threshold = functor->getRefinementThreshold();

	for (ErrorVector::iterator iter = maxErrorSubspaces->begin(); iter != maxErrorSubspaces->end(); ++iter) {
		maxErrorValue = (*iter).getContribPerPoint();
		maxErrorSubspaceIndex = (*iter);

		if (maxErrorValue > functor->start() && fabs(maxErrorValue) >= threshold) {
			createSubspace(storage,maxErrorSubspaceIndex,false);
		}
	}
}

void SubspaceRefinement::freeRefineSubspace(GridStorage* storage, RefinementFunctor* functor) {
	//cout << "in c++ subspace refinement";

	//nothing there => nothing to refine
	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}

	// the functor->getRefinementsNum() largest subspaces should be refined.
	size_t refinements_num = functor->getRefinementsNum();
	//cout << "got refinements num\n";
	//subspaces
	ErrorVector maxSubspaces;

	//accumulate error on refinable subspaces
	//cout << "collecting refinable subspaces\n";
	HashErrorStorage errorStorage(storage->dim());
	collectRefinableSubspaces(storage,functor,&errorStorage);

	//select refinements_num highest indicator subspaces
	selectHighestErrorSubspaces(&errorStorage,refinements_num,&maxSubspaces);

	//refine all subspaces which satisfy the refinement criteria
	//cout << "refining subspaces\n";
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
//		std::cout << "parent level: " << parentLevel << ",parent index: " << parentIndex << std::endl;
		parentSubspace = gridIndex;
		parentSubspace.set(dim,parentLevel,parentIndex);
//		std::cout <<parentSubspace.toString() << std::endl;

		storageIndex = storage->find(&parentSubspace);

		if ( parentLevel>0 && storageIndex == storage->end())
		{
			//parent in dimension dim not found - create the parent subspace
//			cout << "creating parent subspace\n";
			createSubspace(storage,parentSubspace,true);

		}else if (parentLevel>0 && storageIndex->first->isLeaf()) {
			//remove the leaf property from all future parents

			GridPointVector removeLeafProperty;
			index_type newSubspace = storageIndex->first;
			newSubspace.resetIndexVector();
			createAllGridPointsOfSubspace(newSubspace,&removeLeafProperty);

			//cout << "reseting leafs on subspace " << storageIndex->first->toString() << "\n";
			for(GridPointVector::iterator gridPointIter = removeLeafProperty.begin(); gridPointIter != removeLeafProperty.end(); ++gridPointIter)
			{
//				cout << "searching for" << gridPointIter->toString() << "- ";
				storageIndex = storage->find(&(*gridPointIter));
//				cout << (storageIndex == storage->end()) << "; ";
				storageIndex->first->setLeaf(false);
//				cout << storageIndex->first->isLeaf()<< endl;
			}
		}
	}

	//get all indices we want to create.
	GridPointVector gridPointsToCreate;
	createAllGridPointsOfSubspace(gridIndex,&gridPointsToCreate);

//	cout << "creating subspace" << gridIndex.toString() << "\n";

	//create all the grid Points.
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
	//sort by adding indices to a limited size pq.
	std::priority_queue<ErrorType,ErrorVector,smallestErrorFirst> largestErrorObjects;
	size_t size = 0;


	for (HashErrorStorage::grid_map_iterator subspaceErrorIterator = subspaceError->begin(); subspaceErrorIterator != subspaceError->end(); subspaceErrorIterator++)
	{
		if(size < refinements_num)
		{
		  largestErrorObjects.push(*(subspaceErrorIterator->first));
		}

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

	//debug method print errors
	cout << "Selected:\n";
	for (ErrorVector::iterator iter = maxErrorSubspaces->begin(); iter != maxErrorSubspaces->end();++iter) {
		ErrorType error = *iter;
		cout << "Subspace" << error.toString() << std::endl;
	}
}

//void SubspaceRefinement::createGridPointWithParents(GridStorage* storage, index_type& child)
//{
//
//	//search parent in grid storage
//	size_t dim = 0;
//	index_type parent = child;
//
//	//search for parent in each dimension
//	while(dim < storage->dim())
//	{
//		//cout << dim << "<" << storage->dim() << "?";
//		//cout << "- YES\n";
//		//set the level and index vector to the parent's in dimension dim.
//
//		level_t parentLevel;
//		index_t parentIndex;
//		child.getParentLevelAndIndex(&parentLevel, &parentIndex,dim);
//		bool onlevelOne = getParentLevelAndIndex(&parent, dim);
//
//		//cout << "searching for parent: " << parent.toString() << "on level" << dim << "\n";
//
//		//search for parent gridpoint
//		GridStorage::grid_map_iterator parentIter = storage->find(&parent);
//		if(!onlevelOne && !storage->has_key(&child))
//		{
//			if(parentIter != storage->end())
//			{
//				//cout << " - Found. Creating its children\n";
//				//was found. refine.
//				refineGridpoint1D(storage,parent,dim);
//
//				//the point is no longer a leave, since it has at least one child
//				parentIter->first->setLeaf(false);
//
//			}else{
//				//cout << " - NOT Found. Continue Recursion\n";
//				//parent not existing. go search the parent's parent and create it.
//				createGridPointWithParents(storage,parent);
//			}
//		}
//		//parent created in dimension dim. go to next dimension.
//		++dim;
//	}
//}

//void SubspaceRefinement::resetIndexVector(AbstractRefinement::index_type* gridPoint){
//
//	for (size_t dim = 0; dim < gridPoint->dim(); ++dim) {
//		gridPoint->set(dim,gridPoint->getLevel(dim),1);
//	}
//}

void SubspaceRefinement::createAllGridPointsOfSubspace(index_type& subspace,
		GridPointVector* gridPoints) {

	createAllGridPointsOfSubspaceHelper(gridPoints,subspace,0);
}

void SubspaceRefinement::createAllGridPointsOfSubspaceHelper(
		GridPointVector* gridPoints, index_type& storageIndex, size_t dim)
{

	//go through every dimension
	//cout << dim << "<" << storageIndex.dim() << " ?\n";
	if (dim < storageIndex.dim()) {

		//get level of subspace on that dimension
		index_t index;
		level_t level;
		storageIndex.get(dim,level,index);
		index = 1;

		//cout << level << " , " << index << "\n";
		//iterate over all allowed indices on that level, in that dimension
		//cout << index << "<" << static_cast <size_t>( 1 << level) << "?\n";
		while(index < static_cast <size_t>( 1 << level)){

			//set gridpoint's index accordingly
			storageIndex.set(dim,level,index);
			//cout << "currently on  " << dim << ", " << level << ", " << index << "\n";
			//recursive call, so that we iterate over all indices on all levels in all dimensions
			createAllGridPointsOfSubspaceHelper(gridPoints,storageIndex,dim+1);
			//move to next admissible
			index=index+2;
		}

	} else {
		//reached end of recursion scheme. add a new grid point.
//		cout << "vector: adding gridpoint " << storageIndex.toString() << "\n";
		gridPoints->push_back(storageIndex);
	}
}

} /* namespace base */
} /* namespace sg */
