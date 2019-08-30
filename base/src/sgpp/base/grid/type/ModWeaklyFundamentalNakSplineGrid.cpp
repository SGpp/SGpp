// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModWeaklyFundamentalNakSplineGrid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

ModWeaklyFundamentalNakSplineGrid::ModWeaklyFundamentalNakSplineGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SWeaklyFundamentalNakSplineModifiedBase(degree));
}

ModWeaklyFundamentalNakSplineGrid::ModWeaklyFundamentalNakSplineGrid(size_t dim, size_t degree)
    : Grid(dim), generator(storage), degree(degree),
      basis_(new SWeaklyFundamentalNakSplineModifiedBase(degree)) {}

ModWeaklyFundamentalNakSplineGrid::~ModWeaklyFundamentalNakSplineGrid() {}

sgpp::base::GridType ModWeaklyFundamentalNakSplineGrid::getType() {
  return sgpp::base::GridType::ModWeaklyFundamentalNakSpline;
}

SBasis& ModWeaklyFundamentalNakSplineGrid::getBasis() { return *basis_; }

size_t ModWeaklyFundamentalNakSplineGrid::getDegree() { return this->degree; }

Grid* ModWeaklyFundamentalNakSplineGrid::unserialize(std::istream& istr) {
  return new ModWeaklyFundamentalNakSplineGrid(istr);
}

void ModWeaklyFundamentalNakSplineGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModWeaklyFundamentalNakSplineGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
