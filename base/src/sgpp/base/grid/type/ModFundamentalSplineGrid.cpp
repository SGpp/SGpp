// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

ModFundamentalSplineGrid::ModFundamentalSplineGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SFundamentalSplineModifiedBase(degree));
}

ModFundamentalSplineGrid::ModFundamentalSplineGrid(size_t dim, size_t degree)
    : Grid(dim),
      generator(storage),
      degree(degree),
      basis_(new SFundamentalSplineModifiedBase(degree)) {}

ModFundamentalSplineGrid::~ModFundamentalSplineGrid() {}

sgpp::base::GridType ModFundamentalSplineGrid::getType() {
  return sgpp::base::GridType::ModFundamentalSpline;
}

const SBasis& ModFundamentalSplineGrid::getBasis() { return *basis_; }

size_t ModFundamentalSplineGrid::getDegree() { return this->degree; }

Grid* ModFundamentalSplineGrid::unserialize(std::istream& istr) {
  return new ModFundamentalSplineGrid(istr);
}

void ModFundamentalSplineGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModFundamentalSplineGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
