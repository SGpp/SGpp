/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "base/grid/generation/hashmap/HashRefinement.hpp"
#include "base/exception/generation_exception.hpp"

using namespace std;

namespace sg {
namespace base {


void HashRefinement::collectRefinablePoints(GridStorage* storage, RefinementFunctor* functor, size_t refinements_num, size_t* max_indices, RefinementFunctor::value_type* max_values) {
	size_t min_idx = 0;

	// max value equals min value
	RefinementFunctor::value_type max_value = max_values[min_idx];
	//size_t max_index = max_indices[min_idx];

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

			cout << "analyzing point " << index.toString() << "\n";

			// test existence of left child
			index.set(d, source_level + 1, 2 * source_index - 1);
			cout << "testing existence of left child " << index.toString() << "\n";
			child_iter = storage->find(&index);

			// if there no more grid points --> test if we should refine the grid
			if (child_iter == end_iter) {
				RefinementFunctor::value_type current_value = fabs((*functor)(storage, iter->second));

				cout << "is refineable\n";

				if (current_value > max_value) {
					// replace the minimal point in result array, find the new  minimal point
					max_values[min_idx] = current_value;
					max_indices[min_idx] = iter->second;
					min_idx = getIndexOfMin(max_values, refinements_num);
					max_value = max_values[min_idx];
					break;
				}
			}

			// test existance of right child
			index.set(d, source_level + 1, 2 * source_index + 1);
			cout << "testing existence of right child " << index.toString() << "\n";
			child_iter = storage->find(&index);

			if (child_iter == end_iter) {
				RefinementFunctor::value_type current_value = fabs((*functor)(storage, iter->second));
				cout << "is refineable\n";

				if (current_value > max_value) {
					// replace the minimal point in result array, find the new minimal point
					max_values[min_idx] = current_value;
					max_indices[min_idx] = iter->second;
					min_idx = getIndexOfMin(max_values, refinements_num);
					max_value = max_values[min_idx];
					break;
				}
			}

			// reset current grid point in dimension d
			index.set(d, source_level, source_index);
		}
	}
}

//@TODO (lettrich, low) refactor to share code with collect refinable points
void HashRefinement::collectRefinableSubspaces(GridStorage* storage,
											   RefinementFunctor* functor,
											   size_t refinements_num,
											   index_type* maxErrorSubspaces,
											   RefinementFunctor::value_type* maxErrorValues) {

	//cout<< "collecting refinable subspaces\n";
	//storage for accumulated errors on each subspace
	SubspaceError subspaceError;


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

			cout <<"testing existence of left child " << index.toString() <<"\n";

			// if there no more grid points --> test if we should refine the grid
			if (child_iter == end_iter) {
				RefinementFunctor::value_type gridPointError = (*functor)(storage, iter->second);

				//check if already in map
				//cout << "running find operation \n";
				SubspaceError::iterator subSpaceErrorIterator = subspaceError.find(index);

				//not found
				if(subSpaceErrorIterator == subspaceError.end())
				{
					cout << "inserting new subspace " << index.toString() <<"\n";
					//insert new subspace
					index_type indexCopy = index;
					//reset the index vector to one.
					resetIndexVector(&indexCopy);
					//cout << indexCopy.toString() << "\n";
					subspaceError.insert(make_pair(indexCopy, errorContainer(gridPointError)));
				}
				//found
				else
				{

					cout << "adding to old subspace" <<  ((index_type)(subSpaceErrorIterator->first)).toString() << "\n" ;
					//add error to subspace
					//cout << "old " << ((errorContainer) subSpaceErrorIterator->second).toString();
					subSpaceErrorIterator->second += gridPointError;
					//cout << " new " << ((errorContainer) subSpaceErrorIterator->second).toString() << "\n";
				}
				//break;
			}

			// test existance of right child
			index.set(d, source_level + 1, 2 * source_index + 1);
			cout <<"testing existence of right child " << index.toString() <<"\n";
			child_iter = storage->find(&index);

			if (child_iter == end_iter) {
				RefinementFunctor::value_type gridPointError = (*functor)(storage, iter->second);

				//check if already in map
				SubspaceError::iterator subSpaceErrorIterator = subspaceError.find(index);

				//not found
				if(subSpaceErrorIterator == subspaceError.end())
				{
					cout << "inserting new subspace " << index.toString() <<"\n";
					//insert new subspace
					index_type indexCopy = index;
					//reset the index vector to one.
					resetIndexVector(&indexCopy);
					//cout << indexCopy.toString() << "\n";
					subspaceError.insert(make_pair(indexCopy, errorContainer(gridPointError)));
				}
				//found
				else
				{
					cout << "adding to old subspace" <<  ((index_type)(subSpaceErrorIterator->first)).toString() << "\n" ;
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

	cout << "validating adminssebility\n";
	testAdmissibility(&subspaceError);


	cout << "sorting refinable subspaces\n";
	//select the subspaces with highest error


	cout << "Map contains:\n";

	//debug method - print map
	for (SubspaceError::iterator errorIter = subspaceError.begin();errorIter != subspaceError.end(); errorIter++)
	{
		cout << "Subspace" << ((index_type) (errorIter->first)).toString() << " with error " <<  ((errorContainer) errorIter->second).getContribPerPoint() <<"\n";

	}


	//variables for subspaces
	size_t smallestMaxErrorIndex = 0;
	RefinementFunctor::value_type smallestMaxErrorValue = maxErrorValues[smallestMaxErrorIndex];

	//init arrays
	for (size_t i = 0; i < refinements_num; ++i) {
		//maxErrorSubspaces[i] = ;
		maxErrorValues[i] = 0;
	}

	for (SubspaceError::iterator subspaceErrorIterator = subspaceError.begin(); subspaceErrorIterator != subspaceError.end(); subspaceErrorIterator++)
	{
		//run test: when the error is smaller, add. When the errors are equal, test which subspace adds less gridpoints.
		if (((errorContainer) subspaceErrorIterator->second).getContribPerPoint() > smallestMaxErrorValue ||
				(((errorContainer) subspaceErrorIterator->second).getContribPerPoint()== smallestMaxErrorValue &&
					((index_type) (subspaceErrorIterator->first)).getLevelSum() < maxErrorSubspaces[smallestMaxErrorIndex].getLevelSum()
						))
		{
			maxErrorValues[smallestMaxErrorIndex] = ((errorContainer) subspaceErrorIterator->second).getContribPerPoint();
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


void HashRefinement::refineGridpointsCollection(GridStorage* storage, RefinementFunctor* functor, size_t refinements_num, size_t* max_indices, RefinementFunctor::value_type* max_values) {
	RefinementFunctor::value_type max_value;
	size_t max_index;
	// now refine all grid points which satisfy the refinement criteria
	double threshold = functor->getRefinementThreshold();

	for (size_t i = 0; i < refinements_num; i++) {
		max_value = max_values[i];
		max_index = max_indices[i];

		if (max_value > functor->start() && fabs(max_value) >= threshold) {
			refineGridpoint(storage, max_index);
		}
	}
}

void HashRefinement::refineSubspaceCollection(GridStorage* storage, RefinementFunctor* functor, size_t refinements_num, index_type* maxErrorSubspaces, RefinementFunctor::value_type* maxErrorValues) {

	cout << "refining subspace collection\n";
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



void HashRefinement::free_refine(GridStorage* storage, RefinementFunctor* functor) {
	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}

	// the functor->getRefinementsNum() largest grid points should be refined.
	// gather them in an array max_values
	size_t refinements_num = functor->getRefinementsNum();
	// values
	RefinementFunctor::value_type* max_values = new RefinementFunctor::value_type[refinements_num];
	// indices
	size_t* max_indices = new size_t [refinements_num];

	// initialization
	for (size_t i = 0; i < refinements_num; i++) {
		max_values[i] = functor->start();
		max_indices[i] = 0;
	}

	collectRefinablePoints(storage, functor, refinements_num, max_indices, max_values);
	// now refine all grid points which satisfy the refinement criteria
	refineGridpointsCollection(storage, functor, refinements_num, max_indices, max_values);
	delete [] max_values;
	delete [] max_indices;

}


void HashRefinement::freeRefineSubspace(GridStorage* storage, RefinementFunctor* functor) {
	cout << "in c++ subspace refinement";
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


size_t HashRefinement::getNumberOfRefinablePoints(GridStorage* storage) {
	size_t counter = 0;

	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}

	index_type index;
	GridStorage::grid_map_iterator end_iter = storage->end();

	// start iterating over whole grid
	for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
		index = *(iter->first);

		GridStorage::grid_map_iterator child_iter;

		// check for each grid point whether it can be refined (i.e., whether not all children exist yet)
		for (size_t d = 0; d < storage->dim(); d++) {
			index_t source_index;
			level_t source_level;
			index.get(d, source_level, source_index);

			// test existance of left child
			index.set(d, source_level + 1, 2 * source_index - 1);
			child_iter = storage->find(&index);

			// if there no more grid points --> test if we should refine the grid
			if (child_iter == end_iter) {
				counter++;
				break;
			}

			// test existance of right child
			index.set(d, source_level + 1, 2 * source_index + 1);
			child_iter = storage->find(&index);

			if (child_iter == end_iter) {
				counter++;
				break;
			}

			// reset current grid point in dimension d
			index.set(d, source_level, source_index);
		}
	}

	return counter;

}


void HashRefinement::refineGridpoint1D(GridStorage* storage, index_type& index, size_t d) {
	index_t source_index;
	level_t source_level;
	index.get(d, source_level, source_index);
	// generate left child, if necessary
	index.set(d, source_level + 1, 2 * source_index - 1);

	if (!storage->has_key(&index)) {
		index.setLeaf(true);
		createGridpoint(storage, index);
	}

	// generate right child, if necessary
	index.set(d, source_level + 1, 2 * source_index + 1);

	if (!storage->has_key(&index)) {
		index.setLeaf(true);
		createGridpoint(storage, index);
	}

	index.set(d, source_level, source_index);
}


void HashRefinement::refineGridpoint(GridStorage* storage, size_t refine_index) {
	index_type index((*storage)[refine_index]);
	//Sets leaf property of index, which is refined to false
	(*storage)[refine_index]->setLeaf(false);

	// @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
	for (size_t d = 0; d < storage->dim(); d++) {
		refineGridpoint1D(storage, index, d);
	}
}




void HashRefinement::createGridpoint(GridStorage* storage, index_type& index) {
	index_t source_index;
	level_t source_level;

	for (size_t d = 0; d < storage->dim(); d++) {
		createGridpoint1D(index, d, storage, source_index, source_level);
	}

	storage->insert(index);
}

void HashRefinement::createSubspace(GridStorage* storage, index_type& index) {

	//cout <<"creating subspace " << index.toString() << "\n";

	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}
	//start recursive call
	createSubspaceHelper(storage,index,0);
}

//@TODO (lettrich, low) improve performance by iterative algorithm
void HashRefinement::createSubspaceHelper(GridStorage* storage,index_type& storageIndex, size_t dim){

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
		//cout << "creating gridpoint " << storageIndex.toString() << "\n";
		createGridpoint(storage,storageIndex);
	}
}


void HashRefinement::resetIndexVector(AbstractRefinement::index_type* gridPoint){

	for (size_t dim = 0; dim < gridPoint->dim(); ++dim) {
		gridPoint->set(dim,gridPoint->getLevel(dim),1);
	}
}

void HashRefinement::testAdmissibility(SubspaceError* subspaceError) {

}

}
}

