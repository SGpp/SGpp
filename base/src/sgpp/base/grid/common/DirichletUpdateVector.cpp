/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>

namespace sg {

  namespace base {

    DirichletUpdateVector::DirichletUpdateVector(GridStorage* storage):  storage(storage) {
    }

    DirichletUpdateVector::~DirichletUpdateVector() {
    }

    void DirichletUpdateVector::applyDirichletConditions(DataVector& updateVector, DataVector& sourceVector) {
      for (size_t i = 0; i < storage->size(); i++) {
        GridIndex* curPoint = (*storage)[i];

        if (curPoint->isInnerPoint() == false) {
          updateVector.set(i, sourceVector.get(i));
        }
      }
    }

    void DirichletUpdateVector::setBoundariesToZero(DataVector& updateVector) {
      for (size_t i = 0; i < storage->size(); i++) {
        GridIndex* curPoint = (*storage)[i];

        if (curPoint->isInnerPoint() == false) {
          updateVector.set(i, 0.0);
        }
      }
    }

    void DirichletUpdateVector::setInnerPointsToZero(DataVector& updateVector) {
      for (size_t i = 0; i < storage->size(); i++) {
        GridIndex* curPoint = (*storage)[i];

        if (curPoint->isInnerPoint() == true) {
          updateVector.set(i, 0.0);
        }
      }
    }

    void DirichletUpdateVector::multiplyBoundary(DataVector& updateVector, double value) {
      for (size_t i = 0; i < storage->size(); i++) {
        GridIndex* curPoint = (*storage)[i];

        if (curPoint->isInnerPoint() == false) {
          updateVector.set(i, updateVector.get(i)*value);
        }
      }
    }

    void DirichletUpdateVector::multiplyBoundaryVector(DataVector& updateVector, DataVector& factor) {
      for (size_t i = 0; i < storage->size(); i++) {
        GridIndex* curPoint = (*storage)[i];

        if (curPoint->isInnerPoint() == false) {
          updateVector.set(i, updateVector.get(i)* factor.get(i));
        }
      }
    }

    void DirichletUpdateVector::multiply(DataVector& updateVector, double value, bool (*predicate)(GridIndex*, GridStorage*)) {
      for (size_t i = 0; i < storage->size(); i++) {
        GridIndex* curPoint = (*storage)[i];

        if (predicate(curPoint, storage)) {
          updateVector.set(i, updateVector.get(i)*value);
        }
      }
    }

  }
}
