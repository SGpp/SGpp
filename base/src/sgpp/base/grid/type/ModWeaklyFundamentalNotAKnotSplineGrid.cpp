// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModWeaklyFundamentalNotAKnotSplineGrid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

ModWeaklyFundamentalNotAKnotSplineGrid::ModWeaklyFundamentalNotAKnotSplineGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SWeaklyFundamentalNotAKnotSplineModifiedBase(degree));
}

ModWeaklyFundamentalNotAKnotSplineGrid::ModWeaklyFundamentalNotAKnotSplineGrid(size_t dim, size_t degree)
    : Grid(dim), generator(storage), degree(degree),
      basis_(new SWeaklyFundamentalNotAKnotSplineModifiedBase(degree)) {}

ModWeaklyFundamentalNotAKnotSplineGrid::~ModWeaklyFundamentalNotAKnotSplineGrid() {}

sgpp::base::GridType ModWeaklyFundamentalNotAKnotSplineGrid::getType() {
  return sgpp::base::GridType::ModWeaklyFundamentalNotAKnotSpline;
}

SBasis& ModWeaklyFundamentalNotAKnotSplineGrid::getBasis() { return *basis_; }

size_t ModWeaklyFundamentalNotAKnotSplineGrid::getDegree() { return this->degree; }

Grid* ModWeaklyFundamentalNotAKnotSplineGrid::unserialize(std::istream& istr) {
  return new ModWeaklyFundamentalNotAKnotSplineGrid(istr);
}

void ModWeaklyFundamentalNotAKnotSplineGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModWeaklyFundamentalNotAKnotSplineGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
