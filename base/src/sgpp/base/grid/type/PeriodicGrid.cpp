// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/PeriodicGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearPeriodicBasis.hpp>

#include <sgpp/base/grid/generation/PeriodicGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

PeriodicGrid::PeriodicGrid(std::istream& istr) :
  Grid(istr) {
}

PeriodicGrid::PeriodicGrid(size_t dim) :
  Grid(dim) {
}

PeriodicGrid::~PeriodicGrid() {
}

SGPP::base::GridType PeriodicGrid::getType() {
  return SGPP::base::GridType::Periodic;
}

Grid* PeriodicGrid::unserialize(std::istream& istr) {
  return new PeriodicGrid(istr);
}

const SBasis& PeriodicGrid::getBasis() {
  static SLinearPeriodicBasis basis;
  return basis;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* PeriodicGrid::createGridGenerator() {
  return new PeriodicGridGenerator(this->storage);
}


}  // namespace base
}  // namespace SGPP
