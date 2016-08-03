// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

LinearGrid::LinearGrid(std::istream& istr) :
  Grid(istr),
  generator(storage) {
}

LinearGrid::LinearGrid(size_t dim) :
  Grid(dim),
  generator(storage) {
}

LinearGrid::LinearGrid(BoundingBox& BB) :
  Grid(BB),
  generator(storage) {
}

LinearGrid::~LinearGrid() {
}

sgpp::base::GridType LinearGrid::getType() {
  return sgpp::base::GridType::Linear;
}

SBasis& LinearGrid::getBasis() {
  static SLinearBase basis;
  return basis;
}

Grid* LinearGrid::unserialize(std::istream& istr) {
  return new LinearGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& LinearGrid::getGenerator() {
  return generator;
}


}  // namespace base
}  // namespace sgpp
