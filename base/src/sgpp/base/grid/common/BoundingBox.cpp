// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>

#include <sstream>
#include <string>


namespace sgpp {
namespace base {

BoundingBox::BoundingBox(size_t dimension) {
  this->dimension = dimension;
  dimensionBoundaries = new BoundingBox1D[dimension];

  for (size_t i = 0; i < dimension; i++) {
    dimensionBoundaries[i].leftBoundary = 0.0;
    dimensionBoundaries[i].rightBoundary = 1.0;
    dimensionBoundaries[i].bDirichletLeft = false;
    dimensionBoundaries[i].bDirichletRight = false;
  }

  bTrivialCube = true;
}

BoundingBox::BoundingBox(size_t dimension, const BoundingBox1D* boundaries) {
  bTrivialCube = true;
  this->dimension = dimension;
  dimensionBoundaries = new BoundingBox1D[dimension];

  for (size_t i = 0; i < dimension; i++) {
    dimensionBoundaries[i] = boundaries[i];

    if (dimensionBoundaries[i].leftBoundary != 0.0
        || dimensionBoundaries[i].rightBoundary != 1.0) {
      bTrivialCube = false;
    }
  }
}

BoundingBox::BoundingBox(const BoundingBox& copyBoundingBox) {
  bTrivialCube = true;
  dimension = copyBoundingBox.getDimensions();
  dimensionBoundaries = new BoundingBox1D[dimension];

  for (size_t i = 0; i < dimension; i++) {
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

void BoundingBox::setBoundary(size_t d, const BoundingBox1D& newBoundaries) {
  dimensionBoundaries[d] = newBoundaries;

  if (dimensionBoundaries[d].leftBoundary != 0.0
      || dimensionBoundaries[d].rightBoundary != 1.0) {
    bTrivialCube = false;
  }
}

BoundingBox1D BoundingBox::getBoundary(size_t d) const {
  return dimensionBoundaries[d];
}

size_t BoundingBox::getDimensions() const {
  return dimension;
}

double BoundingBox::getIntervalWidth(size_t d) const {
  return dimensionBoundaries[d].rightBoundary - dimensionBoundaries[d].leftBoundary;
}

double BoundingBox::getIntervalOffset(size_t d) const {
  return dimensionBoundaries[d].leftBoundary;
}

bool BoundingBox::isTrivialCube() const {
  return bTrivialCube;
}

bool BoundingBox::hasDirichletBoundaryLeft(size_t d) const {
  return dimensionBoundaries[d].bDirichletLeft;
}

bool BoundingBox::hasDirichletBoundaryRight(size_t d) const {
  return dimensionBoundaries[d].bDirichletRight;
}

void BoundingBox::toString(std::string& text) const {
  std::stringstream str;

  for (size_t d = 0; d < dimension; d++) {
    str << "Dimensions " << d << "(" << dimensionBoundaries[d].leftBoundary
        << "," <<  dimensionBoundaries[d].rightBoundary << ")\n";
  }

  text = str.str();
}

std::string BoundingBox::toString() const {
  std::string str;
  toString(str);
  return str;
}

}  // namespace base
}  // namespace sgpp
