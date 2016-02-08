// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>

#include <sstream>
#include <string>


namespace SGPP {
namespace base {

BoundingBox::BoundingBox(size_t dim) {
  nDim = dim;
  dimensionBoundaries = new DimensionBoundary[nDim];

  for (size_t i = 0; i < nDim; i++) {
    dimensionBoundaries[i].leftBoundary = 0.0;
    dimensionBoundaries[i].rightBoundary = 1.0;
    dimensionBoundaries[i].bDirichletLeft = false;
    dimensionBoundaries[i].bDirichletRight = false;
  }

  bTrivialCube = true;
}

BoundingBox::BoundingBox(size_t dim, DimensionBoundary* boundaries) {
  bTrivialCube = true;
  nDim = dim;
  dimensionBoundaries = new DimensionBoundary[nDim];

  for (size_t i = 0; i < nDim; i++) {
    dimensionBoundaries[i] = boundaries[i];

    if (dimensionBoundaries[i].leftBoundary != 0.0
        || dimensionBoundaries[i].rightBoundary != 1.0) {
      bTrivialCube = false;
    }
  }
}

BoundingBox::BoundingBox(BoundingBox& copyBoundingBox) {
  bTrivialCube = true;
  nDim = copyBoundingBox.getDimensions();
  dimensionBoundaries = new DimensionBoundary[nDim];

  for (size_t i = 0; i < nDim; i++) {
    dimensionBoundaries[i] = copyBoundingBox.getBoundary(i);

    if (dimensionBoundaries[i].leftBoundary != 0.0
        || dimensionBoundaries[i].rightBoundary != 1.0) {
      bTrivialCube = false;
    }
  }
}

BoundingBox::~BoundingBox() {
  delete[] dimensionBoundaries;
}

void BoundingBox::setBoundary(size_t dimension,
                              DimensionBoundary& newBoundaries) {
  dimensionBoundaries[dimension] = newBoundaries;

  if (dimensionBoundaries[dimension].leftBoundary != 0.0
      || dimensionBoundaries[dimension].rightBoundary != 1.0) {
    bTrivialCube = false;
  }
}

DimensionBoundary BoundingBox::getBoundary(size_t dimension) {
  return dimensionBoundaries[dimension];
}

size_t BoundingBox::getDimensions() {
  return nDim;
}

float_t BoundingBox::getIntervalWidth(size_t dimension) {
  return dimensionBoundaries[dimension].rightBoundary -
         dimensionBoundaries[dimension].leftBoundary;
}

float_t BoundingBox::getIntervalOffset(size_t dimension) {
  return dimensionBoundaries[dimension].leftBoundary;
}

bool BoundingBox::isTrivialCube() {
  return bTrivialCube;
}

bool BoundingBox::hasDirichletBoundaryLeft(size_t dimension) {
  return dimensionBoundaries[dimension].bDirichletLeft;
}

bool BoundingBox::hasDirichletBoundaryRight(size_t dimension) {
  return dimensionBoundaries[dimension].bDirichletRight;
}


void BoundingBox::toString(std::string& text) {
  std::stringstream str;

  for (size_t d = 0; d < nDim; d++) {
    str << "Dimensions " << d << "(" << dimensionBoundaries[d].leftBoundary
        << "," <<  dimensionBoundaries[d].rightBoundary << ")\n";
  }

  text = str.str();
}

std::string BoundingBox::toString() {
  std::string str;
  toString(str);
  return str;
}

}  // namespace base
}  // namespace SGPP
