// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/WaveletBoundaryGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

WaveletBoundaryGrid::WaveletBoundaryGrid(std::istream& istr) :
  Grid(istr),
  generator(storage, boundaryLevel),
  boundaryLevel(0) {
  istr >> boundaryLevel;
}

WaveletBoundaryGrid::WaveletBoundaryGrid(size_t dim,
    level_t boundaryLevel) :
  Grid(dim),
  generator(storage, boundaryLevel),
  boundaryLevel(boundaryLevel) {
}

WaveletBoundaryGrid::~WaveletBoundaryGrid() {
}

SGPP::base::GridType WaveletBoundaryGrid::getType() {
  return SGPP::base::GridType::WaveletBoundary;
}

const SBasis& WaveletBoundaryGrid::getBasis() {
  static SWaveletBoundaryBase basis;
  return basis;
}

std::unique_ptr<Grid> WaveletBoundaryGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new WaveletBoundaryGrid(istr));
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& WaveletBoundaryGrid::getGenerator() {
  return generator;
}



}  // namespace base
}  // namespace SGPP
