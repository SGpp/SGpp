// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>

namespace sgpp {
namespace base {

LinearBoundaryGrid::LinearBoundaryGrid(std::istream& istr)
    : Grid(istr), generator(storage, boundaryLevel), boundaryLevel(0) {
  istr >> boundaryLevel;
}

LinearBoundaryGrid::LinearBoundaryGrid(size_t dim, level_t boundaryLevel)
    : Grid(dim), generator(storage, boundaryLevel), boundaryLevel(boundaryLevel) {}

LinearBoundaryGrid::LinearBoundaryGrid(BoundingBox& BB, level_t boundaryLevel)
    : Grid(BB), generator(storage, boundaryLevel), boundaryLevel(boundaryLevel) {}

LinearBoundaryGrid::~LinearBoundaryGrid() {}

sgpp::base::GridType LinearBoundaryGrid::getType() { return sgpp::base::GridType::LinearBoundary; }

const SBasis& LinearBoundaryGrid::getBasis() {
  static SLinearBoundaryBase basis;
  return basis;
}

std::unique_ptr<Grid> LinearBoundaryGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new LinearBoundaryGrid(istr));
}

void LinearBoundaryGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& LinearBoundaryGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
