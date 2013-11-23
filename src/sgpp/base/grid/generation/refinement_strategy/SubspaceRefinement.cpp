/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "SubspaceRefinement.hpp"

using namespace std;


namespace sg {
namespace base {

void SubspaceRefinement::collectRefinableSubspaces(GridStorage* storage,
											   RefinementFunctor* functor,
											   size_t refinements_num,
											   index_type* maxErrorSubspaces,
											   RefinementFunctor::value_type* maxErrorValues) {

	//storage for accumulated errors on each subspace
	SubspaceErrorStorage subspaceError;


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
				RefinementFunctor::value_type gridPointError = (*functor)(storage, iter->second);

				//check if already in map
				//cout << "running find operation \n";
				SubspaceErrorStorage::iterator subSpaceErrorIterator = subspaceError.find(index);

				//not found
				if(subSpaceErrorIterator == subspaceError.end())
				{
					//cout << "inserting new subspace " << index.toString() <<"\n";
					//insert new subspace
					index_type indexCopy = index;
					//reset the index vector to one.
					resetIndexVector(&indexCopy);
					//cout << indexCopy.toString() << "\n";
					subspaceError.insert(make_pair(indexCopy, ErrorContainer(gridPointError)));
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
				RefinementFunctor::value_type gridPointError = (*functor)(storage, iter->second);

				//check if already in map
				SubspaceErrorStorage::iterator subSpaceErrorIterator = subspaceError.find(index);

				//not found
				if(subSpaceErrorIterator == subspaceError.end())
				{
					//cout << "inserting new subspace " << index.toString() <<"\n";
					//insert new subspace
					index_type indexCopy = index;
					//reset the index vector to one.
					resetIndexVector(&indexCopy);
					//cout << indexCopy.toString() << "\n";
					subspaceError.insert(make_pair(indexCopy, ErrorContainer(gridPointError)));
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

	//cout << "validating adminssebility\n";
	//testAdmissibility(&subspaceError);


	//cout << "sorting refinable subspaces\n";
	//select the subspaces with highest error


	cout << "Map contains:\n";

	//debug method - print map
	for (SubspaceErrorStorage::iterator errorIter = subspaceError.begin();errorIter != subspaceError.end(); errorIter++)
	{
		cout << "Subspace" << ((index_type) (errorIter->first)).toString() << " with error " <<  ((ErrorContainer) errorIter->second).getContribPerPoint() <<"\n";

	}


	//variables for subspaces
	size_t smallestMaxErrorIndex = 0;
	RefinementFunctor::value_type smallestMaxErrorValue = maxErrorValues[smallestMaxErrorIndex];

	//init arrays
	for (size_t i = 0; i < refinements_num; ++i) {
		//maxErrorSubspaces[i] = ;
		maxErrorValues[i] = 0;
	}

	for (SubspaceErrorStorage::iterator subspaceErrorIterator = subspaceError.begin(); subspaceErrorIterator != subspaceError.end(); subspaceErrorIterator++)
	{
		//run test: when the error is smaller, add. When the errors are equal, test which subspace adds less gridpoints.
		if (((ErrorContainer) subspaceErrorIterator->second).getContribPerPoint() > smallestMaxErrorValue ||
				(((ErrorContainer) subspaceErrorIterator->second).getContribPerPoint()== smallestMaxErrorValue &&
					((index_type) (subspaceErrorIterator->first)).getLevelSum() < maxErrorSubspaces[smallestMaxErrorIndex].getLevelSum()
						))
		{
			maxErrorValues[smallestMaxErrorIndex] = ((ErrorContainer) subspaceErrorIterator->second).getContribPerPoint();
			maxErrorSubspaces[smallestMaxErrorIndex] = subspaceErrorIterator->first;
			smallestMaxErrorIndex = getIndexOfMin(maxErrorValues, refinements_num);
			smallestMaxErrorValue = maxErrorValues[smallestMaxErrorIndex];
		}
	}


	//debug method print errors
	cout << "Selected:\n";
	for (size_t i = 0; i < refinements_num; ++i) {
		cout << "Subspace" << maxErrorSubspaces[i].toString() << "," << maxErrorValues[i] <<"\n";
	}

}

void SubspaceRefinement::refineSubspaceCollection(GridStorage* storage, RefinementFunctor* functor, size_t refinements_num, index_type* maxErrorSubspaces, RefinementFunctor::value_type* maxErrorValues) {

	//cout << "refining subspace collection\n";
	RefinementFunctor::value_type maxErrorValue;
	index_type maxErrorSubspaceIndex;
	// now refine all grid points which satisfy the refinement criteria
	double threshold = functor->getRefinementThreshold();

	for (size_t i = 0; i < refinements_num; i++) {
		maxErrorValue = maxErrorValues[i];
		maxErrorSubspaceIndex = maxErrorSubspaces[i];

		if (maxErrorValue > functor->start() && fabs(maxErrorValue) >= threshold) {
			createSubspace(storage,maxErrorSubspaceIndex);
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
	cout << "got refinements num\n";
	//subspaces
	index_type* maxSubspaces = new index_type[refinements_num];
	//subspace errors
	RefinementFunctor::value_type* errorsPerSubspace = new RefinementFunctor::value_type[refinements_num];

	//accumulate error on refinable subspaces
	cout << "collecting refinable subspaces\n";
	collectRefinableSubspaces(storage,functor,refinements_num,maxSubspaces,errorsPerSubspace);

	//refine all subspaces which satisfy the refinement criteria
	cout << "refining subspaces\n";
	refineSubspaceCollection(storage,functor,refinements_num,maxSubspaces,errorsPerSubspace);

	delete [] maxSubspaces;
	delete [] errorsPerSubspace;

}

void SubspaceRefinement::createSubspace(GridStorage* storage, index_type& index) {

	cout <<"creating subspace " << index.toString() << "\n";

	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}
	//start recursive call
	createSubspaceHelper(storage,index,0);
}

void SubspaceRefinement::createSubspaceHelper(GridStorage* storage,index_type& storageIndex, size_t dim){

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
			createSubspaceHelper(storage,storageIndex,dim+1);
			//move to next admissible
			index=index+2;
		}

	} else {
		//reached end of recursion scheme. add a new grid point.
		cout << "creating gridpoint " << storageIndex.toString() << "\n";
		//createGridpoint(storage,storageIndex);
		createGridPointWithParents(storage,storageIndex);
	}
}

void SubspaceRefinement::createGridPointWithParents(GridStorage* storage, index_type& child)
{
	//search parent in grid storage
	size_t dim = 0;
	index_type parent = child;

	//search for parent in each dimension
	while(dim < storage->dim())
	{
		//cout << dim << "<" << storage->dim() << "?";
		//cout << "- YES\n";
		//set the level and index vector to the parent's in dimension dim.

		bool onlevelOne = getParentLevelAndIndex(&parent, dim);

		//cout << "searching for parent: " << parent.toString() << "on level" << dim << "\n";

		//search for parent gridpoint
		GridStorage::grid_map_iterator parentIter = storage->find(&parent);
		if(!onlevelOne && !storage->has_key(&child))
		{
			if(parentIter != storage->end())
			{
				//cout << " - Found. Creating its children\n";
				//was found. refine.
				refineGridpoint1D(storage,parent,dim);

				//the point is no longer a leave, since it has at least one child
				parentIter->first->setLeaf(false);

			}else{
				//cout << " - NOT Found. Continue Recursion\n";
				//parent not existing. go search the parent's parent and create it.
				createGridPointWithParents(storage,parent);
			}
		}
		//parent created in dimension dim. go to next dimension.
		++dim;
	}
}

void SubspaceRefinement::resetIndexVector(AbstractRefinement::index_type* gridPoint){

	for (size_t dim = 0; dim < gridPoint->dim(); ++dim) {
		gridPoint->set(dim,gridPoint->getLevel(dim),1);
	}
}

bool SubspaceRefinement::getParentLevelAndIndex(index_type* child,size_t dim)
{
	//get the parents index
	level_t parentLevel = child->getLevel(dim);
	index_t parentIndex = child->getIndex(dim);
	//cout << "analyzing GP " << (*child).toString() << "\n";
	if (parentLevel>1)
	{
		//cout << "parent index = " << parentIndex << ", " << (parentIndex-1)%2 << "\n";
		if(((parentIndex-1)/2)%2 == 0)
		{
			//child is on the left
			child->set(dim,parentLevel-1,(parentIndex+1)/2);
			//cout << "(left) setting child to " << (*child).toString() << "\n";
		}else{
			//child is on the right
			child->set(dim,parentLevel-1,(parentIndex-1)/2);
			//cout << "(right) setting child to " << (*child).toString() << "\n";
		}
		//cout << "child is NOT on level 1\n";
		return false;
	}
	//cout << "child is on level 1\n";
	return true;
}

} /* namespace base */
} /* namespace sg */
