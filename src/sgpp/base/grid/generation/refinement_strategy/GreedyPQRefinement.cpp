/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "GreedyPQRefinement.hpp"
#include "../../../basis/linear/noboundary/LinearBasis.hpp"
#include "../../../basis/linear/boundary/LinearBoundaryBasis.hpp"
#include "../../../basis/modlinear/ModifiedLinearBasis.hpp"

namespace sg {
namespace base {

GreedyPQRefinement::GreedyPQRefinement(AbstractRefinement* refinement,
		DataMatrix* dataSet,
		DataVector* errors,
		Grid* grid): RefinementDecorator(refinement)
{
	this->dataSet=dataSet;
	this->errors=errors;

	basisTypes["linear"] = 0;
	basisTypes["linearBoundary"] = 1;
	basisTypes["modlinear"] = 2;

	//get the name of the grid
	BasisTypes::iterator it = basisTypes.find(grid->getType());
	if (it == basisTypes.end()) {
		//error - gridType not implemented.
		//@TODO: (lettrich), throw exception
	}else {
		gridType = it->second;
	}
}

GreedyPQRefinement::~GreedyPQRefinement() {
}

void GreedyPQRefinement::free_refine(GridStorage* storage,
		RefinementFunctor* functor) {

	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}

	// the functor->getRefinementsNum() largest grid points should be refined.
	// gather them in an array max_values
	size_t refinements_num = functor->getRefinementsNum();
	// values
	RefinementFunctor::value_type* max_values = new RefinementFunctor::value_type[refinements_num];
	// indices
	index_type*  points = new index_type [refinements_num];
	size_t* max_indices = new size_t [refinements_num];

	// initialization
	for (size_t i = 0; i < refinements_num; i++) {
		max_values[i] = functor->start();
	}

	//get all refinable points sorted in multimap.
	collectAllRefinablePointsSorted(storage,functor, &gridObjectsSortedByError);

	//get the points with highest error indicator
	gridObjectsSortedByError.collectHighestErrorIndicators(refinements_num,max_values,points);

	//convert points to grid storage indices
	for (size_t i = 0; i < refinements_num; ++i) {
		max_indices[i]= storage->find(&points[i])->second;
	}

	//refine them
	refineGridpointsCollection(storage,functor,refinements_num,max_indices,max_values);

	//the points were refined. clean up multimap.
	gridObjectsSortedByError.removeAddedErrorIndicators(refinements_num,max_values,points);

	//destroy all created variables.
	delete [] max_values;
	delete [] max_indices;
	delete [] points;
}

void GreedyPQRefinement::freeRefineSubspace(GridStorage* storage,
		RefinementFunctor* functor) {

	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}

	// the functor->getRefinementsNum() largest grid points should be refined.
	// gather them in an array max_values
	size_t refinements_num = functor->getRefinementsNum();
	// values
	RefinementFunctor::value_type* max_values = new RefinementFunctor::value_type[refinements_num];
	// indices
	index_type* max_indices = new index_type [refinements_num];

	// initialization
	for (size_t i = 0; i < refinements_num; i++) {
		max_values[i] = functor->start();
	}
	//collect stuff
	//collectRefinableSubspacesSorted(storage,functor,refinements_num,)

	//get the points with highest error indicator
	gridObjectsSortedByError.collectHighestErrorIndicators(refinements_num,max_values,max_indices);

	//refine them
	refineSubspaceCollection(storage,functor,refinements_num,max_indices, max_values);

	//the points were refined. clean up multimap.
	gridObjectsSortedByError.removeAddedErrorIndicators(refinements_num,max_values,max_indices);

	//destroy all created variables.
	delete [] max_values;
	delete [] max_indices;
}

void GreedyPQRefinement::collectAllRefinablePointsSorted(GridStorage* storage,
		RefinementFunctor* functor,
		SortedGridObjectContainer* sortedRefinablePoints) {


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

			// if there no more grid points --> test if we should refine the grid
			if (child_iter == end_iter) {

				//create a new error container
				GridPointErrorContainer newRefinablePoint;
				//set the parent as error;
				newRefinablePoint.gridPoint = iter->first;

				//this sets the child as the child - easier data structure!?!
				//newRefinablePoint.gridPoint = index;
				//calculate the error indicators
				calculateErrorIndicatorForObject(storage,functor,&newRefinablePoint);
				gridObjectsSortedByError.addErrorObject(newRefinablePoint);
				break;
			}

			// test existance of right child
			index.set(d, source_level + 1, 2 * source_index + 1);
			child_iter = storage->find(&index);

			if (child_iter == end_iter) {

				//create a new error container
				GridPointErrorContainer newRefinablePoint;
				//set the index
				//set the parent as error;
				newRefinablePoint.gridPoint = iter->first;
				//newRefinablePoint.gridPoint = index;
				//calculate the error indicators
				calculateErrorIndicatorForObject(storage,functor,&newRefinablePoint);
				gridObjectsSortedByError.addErrorObject(newRefinablePoint);
				break;
			}

			// reset current grid point in dimension d
			index.set(d, source_level, source_index);
		}
	}
}


//@TODO: (lettrich,low) use functor to calculate the error for different Error Object Types(Grid Points, Subspaces);
void GreedyPQRefinement::calculateErrorIndicatorForObject(
		GridStorage* storage,RefinementFunctor* functor, GridPointErrorContainer* refinableObject)
{
	//get the dimension in which we want to refine
	dimensionType dim = refinableObject->dimension;
	//get the level, and the index on which we are
	level_t level;
	index_t index;
	refinableObject->gridPoint->get(dim,level,index);
	//calculate it's support.
	double supportFloor = (index-1)/(1<<level);
	double supportCeiling= index/(1<<level);
	double supportMid = supportFloor+(supportCeiling-supportFloor)/2;

	//create the indicators for the children
	RefinementFunctor::value_type leftIndicator = 0;
	RefinementFunctor::value_type rightIndicator = 0;

	//helper variables to get basis function value
	vector<double> pointCoordinates;

	//walk through dataSet and select points on support
	for (size_t row = 0; row < dataSet->getNrows(); ++row) {

		//get the values in all dimensions of the current DataPoint
		double valueInDim = dataSet->get(row,dim);

		if (valueInDim>=supportFloor &&
				valueInDim<supportCeiling ) {
			//is on support

			if (valueInDim<supportMid) {
				//belongs to left child
				leftIndicator+=basisFunctionEvalHelper(level+1,2*index-1,valueInDim)*pow(errors->get(row),2);

			}else {
				//belongs to right child
				rightIndicator+=basisFunctionEvalHelper(level+1,2*index+1,valueInDim)*pow(errors->get(row),2);
			}
		}
	}

	//done - write value into the RefinableObject container.
	(*refinableObject)[0]=leftIndicator;
	(*refinableObject)[1]=rightIndicator;
}

double GreedyPQRefinement::basisFunctionEvalHelper(level_t level, index_t index, double value)
{

	switch (gridType) {
		case 0:
			// linear basis
			LinearBasis<level_t,index_t> linBasis;
			return linBasis.eval(level,index,value);
		case 1:
			// linear Basis with Boundaries
			LinearBoundaryBasis<level_t,index_t> linBoundBasis;
			return linBoundBasis.eval(level,index,value);
		case 2:
			// modified linear basis
			ModifiedLinearBasis<level_t,index_t> modLinBasis;
			return modLinBasis.eval(level,index,value);
		default:
			// not found.
			return 0;
	}
}

const SortedGridObjectContainer& GreedyPQRefinement::getGridObjectsSortedByError() const {
	return gridObjectsSortedByError;
}

} /* namespace base */
} /* namespace sg */
