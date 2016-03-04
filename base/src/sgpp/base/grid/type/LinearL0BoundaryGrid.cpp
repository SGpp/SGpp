// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/type/LinearL0BoundaryGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>


namespace sgpp {
namespace base {

LinearL0BoundaryGrid::LinearL0BoundaryGrid(std::istream& istr) :
  Grid(istr),
  generator(storage) {
}

LinearL0BoundaryGrid::LinearL0BoundaryGrid(size_t dim) :
  Grid(dim),
  generator(storage) {
}

LinearL0BoundaryGrid::~LinearL0BoundaryGrid() {
}

sgpp::base::GridType LinearL0BoundaryGrid::getType() {
  return sgpp::base::GridType::LinearL0Boundary;
}

const SBasis& LinearL0BoundaryGrid::getBasis() {
  static SLinearBoundaryBase basis;
  return basis;
}

std::unique_ptr<Grid> LinearL0BoundaryGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new LinearL0BoundaryGrid(istr));
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& LinearL0BoundaryGrid::getGenerator() {
  return generator;
}

}  // namespace base
}  // namespace sgpp
