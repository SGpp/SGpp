/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef EVALCUBOIDGENERATOR_HPP
#define EVALCUBOIDGENERATOR_HPP

#include "data/DataVector.hpp"
#include "grid/common/BoundingBox.hpp"
#include <vector>

namespace sg
{

/**
 * This class builds a cuboid in the d-dimension space. This
 * cuboid is discretesized by a regular full grid.
 *
 * @version $HEAD$
 */
class EvalCuboidGenerator
{
private:
	/// the gris's bounding box
	BoundingBox* myBoundingBox;
	/// number of dimensions
	size_t numDimensions;

	/**
	 * This function is a recursive implementation in order the build the evaluation cuboid
	 *
	 * @param evalPoints vector of dynamic size into which the points are "submitted" during calculation
	 * @param curPoint a current point in the d-dimensional space which which is adjusted during this recursive calculations
	 * @param center the center of the cuboid
	 * @param size the precentage of the whole array the cuboid will cover in a every dimension
	 * @param points number of points used in every dimension
	 */
	void getCuboidEvalPoints(std::vector<DataVector>& evalPoints, DataVector& curPoint, std::vector<double>& center, double size, size_t points, size_t curDim);

public:
	/**
	 * Constructor
	 *
	 * @param BB the grid's bounding box
	 */
	EvalCuboidGenerator(BoundingBox& BB, size_t numDims);

	/**
	 * Destructor
	 */
	~EvalCuboidGenerator();

	/**
	 * This function builds an cuboid which will be stored into the EvaluationPoint
	 * variable of this function.
	 * This is by done by building a cuboid around a given center. The size
	 * of the cuboid is determined in every dimension by a fix percent size of the interval in that dimension.
	 * In addition there is a fix number of EvalutionPoints in every dimension. Be aware that this
	 * function returns point to the power of d points.
	 *
	 * @param EvaluationPoints DataVector that will contain the evaluation points afterwards
	 * @param center the center of the cuboid
	 * @param size the precentage of the whole array the cuboid will cover in a every dimension
	 * @param points number of points used in every dimension
	 */
	void getEvaluationCuboid(DataVector& EvaluationPoints, std::vector<double>& center, double size, size_t points);
};

}

#endif /* EVALCUBOIDGENERATOR_HPP */
