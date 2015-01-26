/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "PredictiveStackANOVARefinement.hpp"

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

void PredictiveStackANOVARefinement::free_refine(GridStorage* storage,RefinementFunctor* functor)
{
	//nothing there => nothing to refine
	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}

	//get the HashErrorStorage from the ErrorStorage DataStructure.
	HashErrorStorage* refinablePointStorage = availableGridPoints.getHashErrorStorage();

	//accumulate error on refinable grid points
	if (firstRefinement) {

		collectRefinablePoints(storage,functor,refinablePointStorage);

		//add all the grid points to a map, sorted by error.
		availableGridPoints.insertAllIntoErrorMap();
		firstRefinement = false;

	}else{
		updateAdmissiblePoints(storage,functor,&addedInLastRefinement,&availableGridPoints);
	}

//	//DEBUG: print all elements
//	std::cout << "\n\navailable gridPoints" << "\n";
//	std::cout << "=================================================" << "\n";
//	for(HashErrorStorage::grid_map_iterator errorIter = refinablePointStorage->begin(); errorIter != refinablePointStorage->end(); ++errorIter)
//	{
//		std::cout << errorIter->first->toString() << std::endl;
//	}


	//-refine all grid points which satisfy the refinement criteria
	//-empty addedInLastRefinement and insert the subspaces that satisfy the refinement criteria into the storage
	//- remove newly created subspaces from available subspaces.
	refineGridpointsCollection(storage,&availableGridPoints,&addedInLastRefinement,functor);


//	//DEBUG: list all points that have remained in map
//	std::cout << "\nremaining gridPoints" << "\n";
//	std::cout << "=================================================" << "\n";
//	for(HashErrorStorage::grid_map_iterator errorIter = refinablePointStorage->begin(); errorIter != refinablePointStorage->end(); ++errorIter)
//	{
//		std::cout << errorIter->first->toString() << std::endl;
//	}
}

void PredictiveStackANOVARefinement::refineGridpointsCollection(GridStorage* storage,
		ErrorStorage* errorStorage,
		ErrorVector* addedInThisStep,
		RefinementFunctor* functor)
{

	addedInThisStep->clear();

	//get the map that contains all the grid points sorted by error indicator.
	ErrorMap* errorMap = errorStorage->getErrorMap();

	// now refine all subspaces which satisfy the refinement criteria
	double threshold = functor->getRefinementThreshold();
	size_t refinements_num = functor->getRefinementsNum();
	//number of refined  grid points in this refinement step
	size_t refined = 0;


	while(refined<refinements_num && !errorMap->empty())
	{
		ErrorType* maxErrorGridPoint = errorStorage->peek();
		//copy of the grid point needed for manipulations
		ErrorType tmp = *maxErrorGridPoint;

		if(maxErrorGridPoint->getContribPerPoint() > functor->start()
				&& fabs(maxErrorGridPoint->getContribPerPoint()) >= threshold)
		{
			createGridpoint(storage,tmp);
			//the parents of the grid point added are no longer leafs
			resetParentLeafs(storage,&tmp);
			++refined;
			//remove the grid point from the list of refinement candidates
			errorStorage->pop();
			//mark the grid point as added in this refinement step.
			addedInThisStep->push_back(*maxErrorGridPoint);
		}

	}

//	//DEBUG: print all grid points that have been added in this refinement step.
//	std::cout << "\nselected grid points" << "\n";
//	std::cout << "=================================================" << "\n";
//	for (ErrorVector::iterator iter = addedInLastStep->begin(); iter != addedInLastStep->end(); ++iter)
//	{
//		std::cout << (*iter).toString() << std::endl;
//	}
}


void PredictiveStackANOVARefinement::updateAdmissiblePoints(
		GridStorage* storage, RefinementFunctor* functor,
		ErrorVector* addedInLastRefinement, ErrorStorage* admissibleGridPoints) {

	//perform dynamic down cast to see, if the functor is a valid predictive refinement indicator;
	PredictiveRefinementIndicator* errorIndicator = dynamic_cast<PredictiveRefinementIndicator*>(functor);

	//get the error storage.
	HashErrorStorage* refinablePoints = admissibleGridPoints->getHashErrorStorage();
	//create a new error storage for all points that are added in this step
	HashErrorStorage newGridPoints(storage->dim());



	//go through all grid points refined in last refinement step.
	for (ErrorVector::iterator addedInLastIter = addedInLastRefinement->begin();
			addedInLastIter != addedInLastRefinement->end(); ++addedInLastIter)
	{
		ErrorType newChild = *addedInLastIter;

		// @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
		// check for each grid point whether it can be refined (i.e., whether not all kids exist yet)
		// if yes, check whether it belongs to the refinements_num largest ones
		for (size_t d = 0; d < storage->dim(); d++) {
			index_t source_index;
			level_t source_level;
			newChild.get(d, source_level, source_index);

			// test existence of left child
			newChild.set(d, source_level + 1, 2 * source_index - 1);
			HashErrorStorage::grid_map_iterator errorIter = refinablePoints->find(&newChild);
			HashErrorStorage::grid_map_iterator newPointsIter = newGridPoints.find(&newChild);

			//the point does not exist, has all parents in all dimensions and is not already amongst the refinement candidates
			if (errorIter == refinablePoints->end()  && newPointsIter == newGridPoints.end() && checkAdmissibility(storage,newChild)) {
				//get indicator value
				RefinementFunctor::value_type errorVal = (*errorIndicator)(&newChild);

				//add the grid point to the refinement candidates, and
				//but do not forget to reset the error values of the previous value!
				newChild.resetError();
				newChild += errorVal;
				newGridPoints.insert(newChild);
			}


			// test existence of right child
			newChild.set(d, source_level + 1, 2 * source_index + 1);
			errorIter = refinablePoints->find(&newChild);
			newPointsIter = newGridPoints.find(&newChild);

			//the point does not exist, has all parents in all dimensions and is not already amongst the refinement candidates.
			if (errorIter == refinablePoints->end()  && newPointsIter == newGridPoints.end() && checkAdmissibility(storage,newChild)) {
				//get indicator value
				RefinementFunctor::value_type errorVal = (*errorIndicator)(&newChild);

				//add the grid point to the refinement candidates
				//but do not forget to reset the error values of the previous value!
				newChild.resetError();
				newChild += errorVal;
				newGridPoints.insert(newChild);
			}

			// reset current grid point in dimension d
			newChild.set(d, source_level, source_index);
		}
	}

	//update available subspaces
	admissibleGridPoints->updateErrors(&newGridPoints);
}



void PredictiveStackANOVARefinement::collectRefinablePoints(
		GridStorage* storage, RefinementFunctor* functor,
		HashErrorStorage* errorStorage) {

	//perform dynamic down cast to see, if the functor is a valid predictive refinement indicator;
	PredictiveRefinementIndicator* errorIndicator = dynamic_cast<PredictiveRefinementIndicator*>(functor);

	//work through all refinable points
	index_type index;
	GridStorage::grid_map_iterator end_iter = storage->end();

	// start iterating over whole grid
	for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
		index = *(iter->first);

		GridStorage::grid_map_iterator child_iter;
		HashErrorStorage::grid_map_iterator errorIter;

		// @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
		// check for each grid point whether it can be refined (i.e., whether not all kids exist yet)
		// if yes, check whether it belongs to the refinements_num largest ones
		for (size_t d = 0; d < storage->dim(); d++) {
			index_t source_index;
			level_t source_level;
			index.get(d, source_level, source_index);

			// test existence of left child and occurrence in error storage
			index.set(d, source_level + 1, 2 * source_index - 1);
			ErrorType errorIndex(index);
			child_iter = storage->find(&index);
			errorIter = errorStorage->find(&errorIndex);

			//the point does not exist, has all parents in all dimensions and is not already amongst the refinement candidates.
			if (child_iter == end_iter && errorIter == errorStorage->end() && checkAdmissibility(storage,index)) {

				//get error for grid point
				RefinementFunctor::value_type errorVal = (*errorIndicator)(&(index));

				//add the grid point to error storage
				errorIndex += errorVal;
				errorStorage->insert(errorIndex);
			}

			// test existence of right child and occurrence in error storage
			index.set(d, source_level + 1, 2 * source_index + 1);
			ErrorType tmp(index);
			errorIndex = tmp;
			child_iter = storage->find(&index);
			errorIter = errorStorage->find(&errorIndex);

			//the point does not exist, has all parents in all dimensions and is not already amongst the refinement candidates.
			if (child_iter == end_iter && errorIter == errorStorage->end() && checkAdmissibility(storage,index)) {

				//get error for grid point
				RefinementFunctor::value_type errorVal = (*errorIndicator)(&(index));

				//add the grid point to error storage
				errorIndex += errorVal;
				errorStorage->insert(errorIndex);
			}

			// reset current grid point in dimension d
			index.set(d, source_level, source_index);
		}
	}

//	//DEBUG - print map
//	std::cout << "ErrorStorage contains:\n";
//	for (HashErrorStorage::grid_map_iterator errorIter = errorStorage->begin();errorIter != errorStorage->end(); errorIter++)
//	{
//		std::cout << "Subspace" << ((ErrorType) (errorIter->first)).toString() << std::endl;
//	}
}

void PredictiveStackANOVARefinement::resetParentLeafs(
		GridStorage* storage, index_type* index) {

	level_t parentLevel;
	index_t parentIndex;
	index_type indexCopy(*index);

	//traverse all parents of a grid point in all dimensions
	for(size_t dim = 0; dim < storage->dim(); ++dim)
	{
		index->getParentLevelAndIndex(&parentLevel,&parentIndex,dim);

		//if parent exists, it is no longer a leaf,
		//and leaf property can be removed
		if(parentLevel>0)
		{
			indexCopy.set(dim,parentLevel,parentIndex);
			indexCopy.setLeaf(false);
			indexCopy = *index;
		}
	}
}

} /* namespace base */
} /* namespace SGPP */

