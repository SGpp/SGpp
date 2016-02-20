// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBoundaryBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearStretchedBoundaryGrid.hpp>
#include <sgpp/base/grid/generation/StretchedBoundaryGridGenerator.hpp>


namespace SGPP {
namespace base {

LinearStretchedBoundaryGrid::LinearStretchedBoundaryGrid(std::istream& istr) :
  Grid(istr) {
}

LinearStretchedBoundaryGrid::LinearStretchedBoundaryGrid(size_t dim) :
  Grid(dim) {
}

LinearStretchedBoundaryGrid::LinearStretchedBoundaryGrid(Stretching& BB) :
  Grid(BB) {
}

LinearStretchedBoundaryGrid::~LinearStretchedBoundaryGrid() {
}

SGPP::base::GridType LinearStretchedBoundaryGrid::getType() {
  return SGPP::base::GridType::LinearStretchedBoundary;
}

const SBasis& LinearStretchedBoundaryGrid::getBasis() {
  static SLinearStretchedBoundaryBase basis;
  return basis;
}

std::unique_ptr<Grid> LinearStretchedBoundaryGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new LinearStretchedBoundaryGrid(istr));
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
std::unique_ptr<GridGenerator> LinearStretchedBoundaryGrid::createGridGenerator() {
  return std::unique_ptr<GridGenerator>(new StretchedBoundaryGridGenerator(*this->storage));
}

}  // namespace base
}  // namespace SGPP
