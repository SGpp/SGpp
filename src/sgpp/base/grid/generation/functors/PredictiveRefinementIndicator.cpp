/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "PredictiveRefinementIndicator.hpp"
#include "base/basis/linear/noboundary/LinearBasis.hpp"
#include "base/basis/linear/boundary/LinearBoundaryBasis.hpp"
#include "base/basis/modlinear/ModifiedLinearBasis.hpp"
#include <map>

namespace sg {
namespace base {


PredictiveRefinementIndicator::PredictiveRefinementIndicator(Grid* grid, DataMatrix* dataSet,DataVector* errorVector,
		size_t refinements_num, double threshold)

{
	//find out what type of grid is used;
	gridType = determineGridType(grid);

	//set global Variables accordingly
	this->dataSet = dataSet;
	this->errorVector = errorVector;
	this->refinementsNum = refinements_num;
	this->threshold = threshold;
}


double PredictiveRefinementIndicator::operator ()(AbstractRefinement::index_type* gridPoint)
{
	double errorIndicator = 0;

		for (size_t dim = 0; dim < gridPoint->dim(); ++dim) {

			AbstractRefinement::level_t level = gridPoint->getLevel(dim);
			AbstractRefinement::index_t index = gridPoint->getIndex(dim);

			double supportFloor = (index-1)/(1<<level);
			double supportCeiling= index/(1<<level);

			//walk through dataSet and select points on support
			for (size_t row = 0; row < dataSet->getNrows(); ++row) {

				//get the values in all dimensions of the current DataPoint
				double valueInDim = dataSet->get(row,dim);

				if (valueInDim>=supportFloor &&
						valueInDim<supportCeiling ) {
					//is on support

					//@TODO:
				errorIndicator += ((RefinementFunctor::value_type) basisFunctionEvalHelper(level,index,valueInDim)*errorVector->get(dim));
				}
			}

		}
		return errorIndicator;

}

double PredictiveRefinementIndicator::operator ()(GridStorage* storage,size_t seq) {
	return errorVector->get(seq);
}


double PredictiveRefinementIndicator::basisFunctionEvalHelper(AbstractRefinement::level_t level, AbstractRefinement::index_t index, double value)
{

	switch (gridType) {
		case 1:
			// linear basis
			LinearBasis<AbstractRefinement::level_t,AbstractRefinement::index_t> linBasis;
			return linBasis.eval(level,index,value);
		case 15:
			// linear Basis with Boundaries
			LinearBoundaryBasis<AbstractRefinement::level_t,AbstractRefinement::index_t> linBoundBasis;
			return linBoundBasis.eval(level,index,value);
		case 8:
			// modified linear basis
			ModifiedLinearBasis<AbstractRefinement::level_t,AbstractRefinement::index_t> modLinBasis;
			return modLinBasis.eval(level,index,value);
		default:
			// not found.
			return 0;
	}
}

size_t PredictiveRefinementIndicator::getRefinementsNum()
{
	return refinementsNum;
}

double PredictiveRefinementIndicator::getRefinementThreshold()
{
	return threshold;
}

double PredictiveRefinementIndicator::start() {
	return 0.0;
}
size_t PredictiveRefinementIndicator::determineGridType(Grid* grid) {

	//define a map where to add the supported grid Types
	typedef std::map<const char*,size_t> GridTypes;
	GridTypes gridTypes;
	// add grid Types and map an integer to them
	gridTypes.insert(std::make_pair("linear",1));
//	gridTypes.insert(std::make_pair("linearstencil",2));
//	gridTypes.insert(std::make_pair("linearStretched",3));
//	gridTypes.insert(std::make_pair("modlinearstencil",4));
//	gridTypes.insert(std::make_pair("modpoly",5));
//	gridTypes.insert(std::make_pair("modWavelet",6));
//	gridTypes.insert(std::make_pair("poly",7));
	gridTypes.insert(std::make_pair("modlinear",8));
//	gridTypes.insert(std::make_pair("modBspline",9));
//	gridTypes.insert(std::make_pair("linearTrapezoidBoundary",10));
//	gridTypes.insert(std::make_pair("linearStretchedTrapezoidBoundary",11));
//	gridTypes.insert(std::make_pair("TruncatedTrapezoid",12));
//	gridTypes.insert(std::make_pair("squareRoot",13));
//	gridTypes.insert(std::make_pair("prewavelet",14));
	gridTypes.insert(std::make_pair("linearBoundary",15));


	//find the integer representation of the grid type and return it.
	//zero otherwise.
	GridTypes::iterator typeIter = gridTypes.find(grid->getType());

	if (typeIter != gridTypes.end()) {
		return typeIter->second;
	}else {
		return 0;
	}

}

} /* namespace base */
} /* namespace sg */
