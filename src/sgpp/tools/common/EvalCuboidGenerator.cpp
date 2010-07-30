/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "tools/common/EvalCuboidGenerator.hpp"

namespace sg
{

EvalCuboidGenerator::EvalCuboidGenerator(BoundingBox& BB, size_t numDims)
{
	myBoundingBox = new BoundingBox(BB);
	numDimensions = numDims;
}

EvalCuboidGenerator::~EvalCuboidGenerator()
{
	delete myBoundingBox;
}

void EvalCuboidGenerator::getCuboidEvalPoints(std::vector<DataVector>& evalPoints, DataVector& curPoint, std::vector<double>& center, double size, size_t points, size_t curDim)
{
	if (curDim == 0)
	{
		if (points > 1)
		{
			for (size_t i = 0; i < points; i++)
			{
				curPoint.set(curDim, std::min(std::max(center[curDim]-(myBoundingBox->getIntervalWidth(curDim)*size)+
						((myBoundingBox->getIntervalWidth(curDim)*size*2/(points-1))*static_cast<double>(i)),
						myBoundingBox->getIntervalOffset(curDim)),
						myBoundingBox->getIntervalOffset(curDim)+myBoundingBox->getIntervalWidth(curDim)));

				evalPoints.push_back(curPoint);
			}
		}
		else
		{
			curPoint.set(curDim, center[curDim]);

			evalPoints.push_back(curPoint);
		}
	}
	else
	{
		if (points > 1)
		{
			for (size_t i = 0; i < points; i++)
			{
				curPoint.set(curDim, std::min(std::max(center[curDim]-(myBoundingBox->getIntervalWidth(curDim)*size)+
						((myBoundingBox->getIntervalWidth(curDim)*size*2/(points-1))*static_cast<double>(i)),
						myBoundingBox->getIntervalOffset(curDim)),
						myBoundingBox->getIntervalOffset(curDim)+myBoundingBox->getIntervalWidth(curDim)));

				getCuboidEvalPoints(evalPoints, curPoint, center, size, points, curDim-1);
			}
		}
		else
		{
			curPoint.set(curDim, center[curDim]);

			getCuboidEvalPoints(evalPoints, curPoint, center, size, points, curDim-1);
		}
	}
}

void EvalCuboidGenerator::getEvaluationCuboid(DataVector& EvaluationPoints, std::vector<double>& center, double size, size_t points)
{
	std::vector<DataVector> evalPoints;
	DataVector curPoint(numDimensions);

	getCuboidEvalPoints(evalPoints, curPoint, center, size, points, numDimensions-1);

	size_t numEvalPoints = evalPoints.size();
	EvaluationPoints.resize(numEvalPoints);

	for (size_t i = 0; i < numEvalPoints; i++)
	{
		EvaluationPoints.setRow(i, evalPoints[i]);
	}
}

}
