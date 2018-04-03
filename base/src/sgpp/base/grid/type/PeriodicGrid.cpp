// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/PeriodicGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearPeriodicBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

PeriodicGrid::PeriodicGrid(std::istream& istr) :
  Grid(istr),
  generator(storage) {
}

PeriodicGrid::PeriodicGrid(size_t dim) :
  Grid(dim),
  generator(storage) {
}

PeriodicGrid::~PeriodicGrid() {
}

sgpp::base::GridType PeriodicGrid::getType() {
  return sgpp::base::GridType::Periodic;
}

Grid* PeriodicGrid::unserialize(std::istream& istr) {
  return new PeriodicGrid(istr);
}

SBasis& PeriodicGrid::getBasis() {
  static SLinearPeriodicBasis basis;
  return basis;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& PeriodicGrid::getGenerator() {
  return generator;
}


}  // namespace base
}  // namespace sgpp
