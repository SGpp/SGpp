/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "base/tools/EvalCuboidGeneratorForStretching.hpp"

namespace sg {
  namespace base {

    EvalCuboidGeneratorForStretching::EvalCuboidGeneratorForStretching() {
    }

    EvalCuboidGeneratorForStretching::~EvalCuboidGeneratorForStretching() {
    }

    void EvalCuboidGeneratorForStretching::getCuboidEvalPoints(std::vector<DataVector>& evalPoints, DataVector& curPoint, Stretching& myStretching, size_t points, size_t curDim) {
      if (curDim == 0) {
        if (points > 1) {
          for (size_t i = 0; i < points; i++) {
            double inc = (myStretching.getBoundary(curDim).rightBoundary - myStretching.getBoundary(curDim).leftBoundary) / static_cast<double>(points - 1);

            curPoint.set(curDim, myStretching.getBoundary(curDim).leftBoundary + (inc * static_cast<double>(i)));

            evalPoints.push_back(curPoint);
          }
        } else {
          curPoint.set(curDim, (myStretching.getBoundary(curDim).leftBoundary + myStretching.getBoundary(curDim).rightBoundary) / 2.0 );

          evalPoints.push_back(curPoint);
        }
      } else {
        if (points > 1) {
          for (size_t i = 0; i < points; i++) {
            double inc = (myStretching.getBoundary(curDim).rightBoundary - myStretching.getBoundary(curDim).leftBoundary) / static_cast<double>(points - 1);

            curPoint.set(curDim, myStretching.getBoundary(curDim).leftBoundary + (inc * static_cast<double>(i)));

            getCuboidEvalPoints(evalPoints, curPoint, myStretching, points, curDim - 1);
          }
        } else {
          curPoint.set(curDim, (myStretching.getBoundary(curDim).leftBoundary + myStretching.getBoundary(curDim).rightBoundary) / 2.0 );

          getCuboidEvalPoints(evalPoints, curPoint, myStretching, points, curDim - 1);
        }
      }
    }

    void EvalCuboidGeneratorForStretching::getEvaluationCuboid(DataMatrix& EvaluationPoints, Stretching& SubDomain, size_t points) {
      std::vector<DataVector> evalPoints;
      DataVector curPoint(SubDomain.getDimensions());

      getCuboidEvalPoints(evalPoints, curPoint, SubDomain, points, SubDomain.getDimensions() - 1);

      size_t numEvalPoints = evalPoints.size();
      EvaluationPoints.resize(numEvalPoints);

      for (size_t i = 0; i < numEvalPoints; i++) {
        EvaluationPoints.setRow(i, evalPoints[i]);
      }
    }

  }
}
