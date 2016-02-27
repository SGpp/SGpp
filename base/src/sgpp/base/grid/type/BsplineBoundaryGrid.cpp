// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/BsplineBoundaryGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

BsplineBoundaryGrid::BsplineBoundaryGrid(std::istream& istr) :
  Grid(istr),
  generator(storage, boundaryLevel),
  degree(1 << 16),
  boundaryLevel(0) {
  istr >> degree;
  istr >> boundaryLevel;
  basis_.reset(new SBsplineBoundaryBase(degree));
}


BsplineBoundaryGrid::BsplineBoundaryGrid(size_t dim,
    size_t degree,
    level_t boundaryLevel) :
  Grid(dim),
  generator(storage, boundaryLevel),
  degree(degree),
  basis_(new SBsplineBoundaryBase(degree)),
  boundaryLevel(boundaryLevel) {
}

BsplineBoundaryGrid::~BsplineBoundaryGrid() {
}

sgpp::base::GridType BsplineBoundaryGrid::getType() {
  return sgpp::base::GridType::BsplineBoundary;
}

const SBasis& BsplineBoundaryGrid::getBasis() {
  return *basis_;
}

size_t BsplineBoundaryGrid::getDegree() {
  return this->degree;
}

std::unique_ptr<Grid> BsplineBoundaryGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new BsplineBoundaryGrid(istr));
}

void BsplineBoundaryGrid::serialize(std::ostream& ostr) {
  this->Grid::serialize(ostr);
  ostr << degree << std::endl;
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& BsplineBoundaryGrid::getGenerator() {
  return generator;
}

}  // namespace base
}  // namespace sgpp
