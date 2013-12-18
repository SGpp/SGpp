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
	//std::cout << "gridType = " << grid->getType() << " corresponding to " << gridType << "\n";

	//set global Variables accordingly
	this->dataSet = dataSet;
	this->errorVector = errorVector;
	this->refinementsNum = refinements_num;
	this->threshold = threshold;
}

double PredictiveRefinementIndicator::operator ()(AbstractRefinement::index_type* gridPoint)
{
	//calculate the floor and ceiling of the support on dimension of the grid point.
	DataVector floorMask(dataSet->getNcols());
	DataVector ceilingMask(dataSet->getNcols());
	buildGPSupportMask(gridPoint,&floorMask,&ceilingMask);

	//level and index vectors for grid point needed below
	AbstractRefinement::level_t level = 0;
	AbstractRefinement::index_t index = 0;
	//value of the data set in the dimension DIM (needed below)
	double valueInDim;

	double errorIndicator = 0;
	//DEBUG
	size_t counter = 0;

	//go through the whole dataset.
	for(size_t row = 0; row < dataSet->getNrows(); ++row)
	{
		//if on the support of the grid point in all dim
		if(isOnSupport(&floorMask,&ceilingMask,row))
		{
			//calculate error Indicator
			for(size_t dim = 0; dim < gridPoint->dim(); ++dim)
			{
				++counter;
				level = gridPoint->getLevel(dim);
				index = gridPoint->getIndex(dim);

				valueInDim = dataSet->get(row,dim);
				errorIndicator += ((RefinementFunctor::value_type) basisFunctionEvalHelper(level,index,valueInDim)*errorVector->get(row));
			}
		}
	}
	std::cout << gridPoint->toString() << " with error estimate " << errorIndicator << ",caused by " << counter << "contribs - in average: " << (errorIndicator/static_cast<double>(counter)) << "\n";
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
	typedef std::map<std::string,size_t> GridTypes;
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
	GridTypes::iterator typeIter = gridTypes.find(std::string(grid->getType()));

	if (typeIter != gridTypes.end()) {
		return typeIter->second;
	}else {
		return 0;
	}
}

bool PredictiveRefinementIndicator::isOnSupport(
		DataVector* floorMask, DataVector* ceilingMask, size_t entry)
{
	for(size_t col = 0; col < dataSet->getNcols(); ++col)
	{
		double valueInDim = dataSet->get(entry,col);
		if(valueInDim < floorMask->get(col) || valueInDim >= ceilingMask->get(col) )
		{
			return false;
		}
	}
	return true;
}

void PredictiveRefinementIndicator::buildGPSupportMask(
		AbstractRefinement::index_type* gridPoint, DataVector* floorMask, DataVector* ceilingMask)
{

	AbstractRefinement::level_t level;
	AbstractRefinement::index_t index;

	for(size_t dim=0; dim<gridPoint->dim();++dim)
	{
		level = gridPoint->getLevel(dim);
		index = gridPoint->getIndex(dim);

		floorMask->set(dim,(index-1.0)/(1<<level));
		ceilingMask->set(dim,index*1.0/(1<<level));
	}
}




} /* namespace base */
} /* namespace sg */
