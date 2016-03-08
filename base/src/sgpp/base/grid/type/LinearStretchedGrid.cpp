// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearStretchedGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

LinearStretchedGrid::LinearStretchedGrid(std::istream& istr) :
  Grid(istr),
  generator(storage) {
}

LinearStretchedGrid::LinearStretchedGrid(size_t dim) :
  Grid(dim),
  generator(storage) {
}

LinearStretchedGrid::LinearStretchedGrid(Stretching& BB) :
  Grid(BB),
  generator(storage) {
}

LinearStretchedGrid::~LinearStretchedGrid() {
}

sgpp::base::GridType LinearStretchedGrid::getType() {
  return sgpp::base::GridType::LinearStretched;
}

const SBasis& LinearStretchedGrid::getBasis() {
  static SLinearStretchedBase basis;
  return basis;
}

std::unique_ptr<Grid> LinearStretchedGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new LinearStretchedGrid(istr));
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& LinearStretchedGrid::getGenerator() {
  return generator;
}

}  // namespace base
}  // namespace sgpp
