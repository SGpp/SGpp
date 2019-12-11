// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/KinkLinearGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearKinkedBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

KinkLinearGrid::KinkLinearGrid(std::istream& istr) :
  Grid(istr),
  generator(storage) {
}

KinkLinearGrid::KinkLinearGrid(size_t dim) :
  Grid(dim),
  generator(storage) {
}

KinkLinearGrid::~KinkLinearGrid() {
}

sgpp::base::GridType KinkLinearGrid::getType() {
  return sgpp::base::GridType::KinkLinear;
}

SBasis& KinkLinearGrid::getBasis() {
  static SLinearKinkedBase basis;
  return basis;
}

Grid* KinkLinearGrid::unserialize(std::istream& istr) {
  return new KinkLinearGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& KinkLinearGrid::getGenerator() {
  return generator;
}


}  // namespace base
}  // namespace sgpp
