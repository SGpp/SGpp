// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>

#include <sstream>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

BoundingBox::BoundingBox(size_t dimension) :
    dimension(dimension),
    boundingBox1Ds(dimension, {0.0, 1.0}) {
}

BoundingBox::BoundingBox(const std::vector<BoundingBox1D>& boundingBox1Ds) :
    dimension(boundingBox1Ds.size()),
    boundingBox1Ds(boundingBox1Ds) {
}

BoundingBox::~BoundingBox() {
}

void BoundingBox::setBoundary(size_t d, const BoundingBox1D& boundingBox1D) {
  boundingBox1Ds[d] = boundingBox1D;
}

size_t BoundingBox::getDimensions() const {
  return dimension;
}

bool BoundingBox::isUnitCube() const {
  for (size_t d = 0; d < dimension; d++) {
    if ((boundingBox1Ds[d].leftBoundary != 0.0) ||
        (boundingBox1Ds[d].rightBoundary != 1.0)) {
      return false;
    }
  }

  return true;
}

bool BoundingBox::hasDirichletBoundaryLeft(size_t d) const {
  return boundingBox1Ds[d].bDirichletLeft;
}

bool BoundingBox::hasDirichletBoundaryRight(size_t d) const {
  return boundingBox1Ds[d].bDirichletRight;
}

std::string BoundingBox::serialize(int version) const {
  std::ostringstream ostream;
  serialize(ostream, version);
  return ostream.str();
}

void BoundingBox::serialize(std::ostream& ostream, int version) const {
  for (size_t d = 0; d < dimension; d++) {
    const BoundingBox1D& bb1D = boundingBox1Ds[d];
    ostream << std::scientific << bb1D.leftBoundary << " " << bb1D.rightBoundary << " "
            << bb1D.bDirichletLeft << " " << bb1D.bDirichletRight << " ";
  }

  ostream << std::endl;
}

void BoundingBox::unserialize(const std::string& istr, int version) {
  std::istringstream istream;
  istream.str(istr);
  unserialize(istream, version);
}

void BoundingBox::unserialize(std::istream& istr, int version) {
  for (size_t d = 0; d < dimension; d++) {
    BoundingBox1D& bb1D = boundingBox1Ds[d];

    istr >> bb1D.leftBoundary;
    istr >> bb1D.rightBoundary;
    istr >> bb1D.bDirichletLeft;
    istr >> bb1D.bDirichletRight;
  }
}

void BoundingBox::toString(std::string& text) const {
  std::stringstream str;

  for (size_t d = 0; d < dimension; d++) {
    str << "Dimensions " << d << "(" << boundingBox1Ds[d].leftBoundary
        << "," <<  boundingBox1Ds[d].rightBoundary << ")\n";
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
