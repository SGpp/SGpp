// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

namespace base {

DirichletUpdateVector::DirichletUpdateVector(GridStorage* storage):  storage(
    storage) {
}

DirichletUpdateVector::~DirichletUpdateVector() {
}

void DirichletUpdateVector::applyDirichletConditions(DataVector& updateVector,
    DataVector& sourceVector) {
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

void DirichletUpdateVector::multiplyBoundary(DataVector& updateVector,
    float_t value) {
  for (size_t i = 0; i < storage->size(); i++) {
    GridIndex* curPoint = (*storage)[i];

    if (curPoint->isInnerPoint() == false) {
      updateVector.set(i, updateVector.get(i)*value);
    }
  }
}

void DirichletUpdateVector::multiplyBoundaryVector(DataVector& updateVector,
    DataVector& factor) {
  for (size_t i = 0; i < storage->size(); i++) {
    GridIndex* curPoint = (*storage)[i];

    if (curPoint->isInnerPoint() == false) {
      updateVector.set(i, updateVector.get(i)* factor.get(i));
    }
  }
}

void DirichletUpdateVector::multiply(DataVector& updateVector, float_t value,
                                     bool (*predicate)(GridIndex*, GridStorage*)) {
  for (size_t i = 0; i < storage->size(); i++) {
    GridIndex* curPoint = (*storage)[i];

    if (predicate(curPoint, storage)) {
      updateVector.set(i, updateVector.get(i)*value);
    }
  }
}

}
}