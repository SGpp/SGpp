// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>


#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearTruncatedBoundaryGrid.hpp>
#include <sgpp/base/grid/generation/GeneralizedBoundaryGridGenerator.hpp>


namespace SGPP {
namespace base {

LinearTruncatedBoundaryGrid::LinearTruncatedBoundaryGrid(std::istream& istr) :
  Grid(istr) {
}

LinearTruncatedBoundaryGrid::LinearTruncatedBoundaryGrid(size_t dim) :
  Grid(dim) {
}

LinearTruncatedBoundaryGrid::LinearTruncatedBoundaryGrid(BoundingBox& BB) :
  Grid(BB) {
}

LinearTruncatedBoundaryGrid::~LinearTruncatedBoundaryGrid() {
}

SGPP::base::GridType LinearTruncatedBoundaryGrid::getType() {
  return SGPP::base::GridType::LinearTruncatedBoundary;
}

const SBasis& LinearTruncatedBoundaryGrid::getBasis() {
  static SLinearBoundaryBase basis;
  return basis;
}

std::unique_ptr<Grid> LinearTruncatedBoundaryGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new LinearTruncatedBoundaryGrid(istr));
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
std::unique_ptr<GridGenerator> LinearTruncatedBoundaryGrid::createGridGenerator() {
  return std::unique_ptr<GridGenerator>(new GeneralizedBoundaryGridGenerator(*this->storage));
}

}  // namespace base
}  // namespace SGPP
