// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearTruncatedBoundaryGrid.hpp>


namespace sgpp {
namespace base {

LinearTruncatedBoundaryGrid::LinearTruncatedBoundaryGrid(std::istream& istr) :
  Grid(istr),
  generator(storage) {
}

LinearTruncatedBoundaryGrid::LinearTruncatedBoundaryGrid(size_t dim) :
  Grid(dim),
  generator(storage) {
}

LinearTruncatedBoundaryGrid::LinearTruncatedBoundaryGrid(BoundingBox& BB) :
  Grid(BB),
  generator(storage) {
}

LinearTruncatedBoundaryGrid::~LinearTruncatedBoundaryGrid() {
}

sgpp::base::GridType LinearTruncatedBoundaryGrid::getType() {
  return sgpp::base::GridType::LinearTruncatedBoundary;
}

const SBasis& LinearTruncatedBoundaryGrid::getBasis() {
  static SLinearBoundaryBase basis;
  return basis;
}

Grid* LinearTruncatedBoundaryGrid::unserialize(std::istream& istr) {
  return new LinearTruncatedBoundaryGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& LinearTruncatedBoundaryGrid::getGenerator() {
  return generator;
}

}  // namespace base
}  // namespace sgpp
