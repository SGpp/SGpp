/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef PREDICTIVEREFINEMENTINDICATOR_HPP_
#define PREDICTIVEREFINEMENTINDICATOR_HPP_

#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "RefinementFunctor.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/datatypes/DataVector.hpp"


namespace sg {
namespace base {

class PredictiveRefinementIndicator: public RefinementFunctor {
public:
	PredictiveRefinementIndicator(Grid* grid, DataMatrix* dataSet,DataVector* errorVector,
			size_t refinements_num = 1, double threshold = 0.0);


	double operator()(AbstractRefinement::index_type* gridPoint);

	virtual double operator()(GridStorage* storage, size_t seq);

	virtual size_t getRefinementsNum();

	virtual double getRefinementThreshold();

	virtual double start();

protected:

	// for each Point in the dataSet, this Contains the squared absolute offset between sparse grid and data point value;
	DataVector* errorVector;

	/// number of grid points to refine
	size_t refinementsNum;
	/// threshold, only the points with greater to equal absolute values of the refinement criterion (e.g. alpha or error) will be refined
	double threshold;

	// data set that will be evaluated
	DataMatrix* dataSet;

	double basisFunctionEvalHelper(AbstractRefinement::level_t level, AbstractRefinement::index_t index, double value);

private:

	size_t determineGridType(Grid* grid);

	size_t gridType;
};

} /* namespace base */
} /* namespace sg */
#endif /* PREDICTIVEREFINEMENTINDICATOR_HPP_ */
