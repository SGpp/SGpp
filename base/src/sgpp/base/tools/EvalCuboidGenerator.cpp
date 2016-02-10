// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/EvalCuboidGenerator.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace SGPP {
namespace base {

EvalCuboidGenerator::EvalCuboidGenerator() {
}

EvalCuboidGenerator::~EvalCuboidGenerator() {
}

void EvalCuboidGenerator::getCuboidEvalPoints(std::vector<DataVector>&
    evalPoints, DataVector& curPoint, BoundingBox& myBoundingBox, size_t points,
    size_t curDim) {
  if (curDim == 0) {
    if (points > 1) {
      for (size_t i = 0; i < points; i++) {
        float_t inc = (myBoundingBox.getBoundary(curDim).rightBoundary -
                       myBoundingBox.getBoundary(curDim).leftBoundary) /
                      static_cast<float_t>(points - 1);

        curPoint.set(curDim, myBoundingBox.getBoundary(curDim).leftBoundary +
                     (inc * static_cast<float_t>(i)));

        evalPoints.push_back(curPoint);
      }
    } else {
      curPoint.set(curDim, (myBoundingBox.getBoundary(curDim).leftBoundary +
                            myBoundingBox.getBoundary(curDim).rightBoundary) /
                   2.0);

      evalPoints.push_back(curPoint);
    }
  } else {
    if (points > 1) {
      for (size_t i = 0; i < points; i++) {
        float_t inc = (myBoundingBox.getBoundary(curDim).rightBoundary -
                       myBoundingBox.getBoundary(curDim).leftBoundary) /
                      static_cast<float_t>(points - 1);

        curPoint.set(curDim, myBoundingBox.getBoundary(curDim).leftBoundary +
                     (inc * static_cast<float_t>(i)));

        getCuboidEvalPoints(evalPoints, curPoint, myBoundingBox, points,
                            curDim - 1);
      }
    } else {
      curPoint.set(curDim, (myBoundingBox.getBoundary(curDim).leftBoundary +
                            myBoundingBox.getBoundary(curDim).rightBoundary) /
                   2.0);

      getCuboidEvalPoints(evalPoints, curPoint, myBoundingBox, points,
                          curDim - 1);
    }
  }
}

void EvalCuboidGenerator::getEvaluationCuboid(DataMatrix& EvaluationPoints,
    BoundingBox& SubDomain, size_t points) {
  std::vector<DataVector> evalPoints;
  DataVector curPoint(SubDomain.getDimensions());

  getCuboidEvalPoints(evalPoints, curPoint, SubDomain, points,
                      SubDomain.getDimensions() - 1);

  size_t numEvalPoints = evalPoints.size();
  EvaluationPoints.resize(numEvalPoints);

  for (size_t i = 0; i < numEvalPoints; i++) {
    EvaluationPoints.setRow(i, evalPoints[i]);
  }
}

}  // namespace base
}  // namespace SGPP
