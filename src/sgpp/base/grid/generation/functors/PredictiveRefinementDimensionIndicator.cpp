/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "PredictiveRefinementDimensionIndicator.hpp"
#include "base/basis/linear/noboundary/LinearBasis.hpp"
#include "base/basis/linear/boundary/LinearBoundaryBasis.hpp"
#include "base/basis/modlinear/ModifiedLinearBasis.hpp"
#include <map>
#include <string>
#include <utility>

namespace sg {
namespace base {


PredictiveRefinementDimensionIndicator::PredictiveRefinementDimensionIndicator(Grid* grid, DataMatrix* dataSet,DataVector* errorVector,
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

double PredictiveRefinementDimensionIndicator::operator ()(AbstractRefinement::index_type* gridPoint)
{
  //TODO make it configurable
  unsigned short int MIN_POINTS_SUPPORT = 0;
  //calculate the floor and ceiling of the support on dimension of the grid point.#
	DataVector floorMask(dataSet->getNcols());
	DataVector ceilingMask(dataSet->getNcols());
	buildGPSupportMask(gridPoint,&floorMask,&ceilingMask);

	//level, index and evaulation of a gridPoint in dimension d
	AbstractRefinement::level_t level = 0;
	AbstractRefinement::index_t index = 0;
	double valueInDim;

	//the actuall value of the errorIndicator
	double errorIndicator = 0.0;
	double denominator = 0.0;
	double r22 = 0.0;
	double r2phi = 0.0;


	//counter of contributions - for DEBUG purposes
	size_t counter = 0;

	//go through the whole dataset. -> if data point on the support of the grid point in all dim then calculate error Indicator.
	for(size_t row = 0; row < dataSet->getNrows(); ++row)
	{
		//if on the support of the grid point in all dim
		if(isOnSupport(&floorMask,&ceilingMask,row))
		{
			//counter for DEBUG
			//++counter;*****
		  double funcval = 1.0;
			//calculate error Indicator
			for(size_t dim = 0; dim < gridPoint->dim(); ++dim)
			{

				level = gridPoint->getLevel(dim);
				index = gridPoint->getIndex(dim);

				valueInDim = dataSet->get(row,dim);

				funcval *=  basisFunctionEvalHelper(level,index,valueInDim);
			}
		  errorIndicator += funcval*errorVector->get(row)/**errorVector->get(row)*/;
		  r22 += errorVector->get(row)*errorVector->get(row);
		  r2phi += funcval * errorVector->get(row);
		  denominator += funcval * funcval;
		  counter++;
		}
	}

	if (denominator != 0 && counter >= MIN_POINTS_SUPPORT){
	  // to match with OnlineRefDim, use this:
	  return (errorIndicator * errorIndicator) / denominator;

	  double a = (errorIndicator/denominator);
	  return /*2*r22 - 2*a*r2phi + a*a*denominator*/ a*(2*r2phi - a*denominator);
	  //return fabs(a);
	}
	else {
	  return 0.0;
	}

}


double PredictiveRefinementDimensionIndicator::operator ()(GridStorage* storage,size_t seq) {
  return errorVector->get(seq);
}


double PredictiveRefinementDimensionIndicator::runOperator(GridStorage* storage,size_t seq) {
	return (*this)(storage->get(seq));
}


double PredictiveRefinementDimensionIndicator::basisFunctionEvalHelper(AbstractRefinement::level_t level, AbstractRefinement::index_t index, double value)
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

size_t PredictiveRefinementDimensionIndicator::getRefinementsNum()
{
	return refinementsNum;
}

double PredictiveRefinementDimensionIndicator::getRefinementThreshold()
{
	return threshold;
}

double PredictiveRefinementDimensionIndicator::start() {
	return 0.0;
}

size_t PredictiveRefinementDimensionIndicator::determineGridType(Grid* grid) {

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

bool PredictiveRefinementDimensionIndicator::isOnSupport(
		DataVector* floorMask, DataVector* ceilingMask, size_t row)
{

	//go through all cols of the dataset
	//=> go through all samples in dataset and check if in dim "col" "valueInDim" is on support
	for(size_t col = 0; col < dataSet->getNcols(); ++col)
	{
		double valueInDim = dataSet->get(row,col);
		if(valueInDim < floorMask->get(col) || valueInDim >= ceilingMask->get(col) )
		{
			return false;
		}
	}
	DataVector vector(dataSet->getNcols());
	dataSet->getRow(row,vector);

	//Debug
//	std::cout << vector.toString() << " is on support of " << floorMask->toString() << " & " << ceilingMask->toString() << std::endl;

	return true;
}

void PredictiveRefinementDimensionIndicator::buildGPSupportMask(
		AbstractRefinement::index_type* gridPoint, DataVector* floorMask, DataVector* ceilingMask)
{

	AbstractRefinement::level_t level;
	AbstractRefinement::index_t index;

	//in each dimension, get level and index, calculate min and max of supp(GridPointBasisFunction)
	for(size_t dim=0; dim<gridPoint->dim();++dim)
	{
		level = gridPoint->getLevel(dim);
		index = gridPoint->getIndex(dim);

		floorMask->set(dim,(index-1.0)/(1<<(level)));
		ceilingMask->set(dim,(index+1.0)/(1<<(level)));

		//DEBUG
//		std::cout << "floor: " << floorMask->get(dim) << "ceiling " << ceilingMask->get(dim) <<std::endl;
	}
}




} /* namespace base */
} /* namespace sg */
