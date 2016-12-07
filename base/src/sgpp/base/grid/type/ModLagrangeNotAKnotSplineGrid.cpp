// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModLagrangeNotAKnotSplineGrid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

ModLagrangeNotAKnotSplineGrid::ModLagrangeNotAKnotSplineGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SLagrangeNotAKnotSplineModifiedBase(degree));
}

ModLagrangeNotAKnotSplineGrid::ModLagrangeNotAKnotSplineGrid(size_t dim, size_t degree)
    : Grid(dim), generator(storage), degree(degree),
      basis_(new SLagrangeNotAKnotSplineModifiedBase(degree)) {}

ModLagrangeNotAKnotSplineGrid::~ModLagrangeNotAKnotSplineGrid() {}

sgpp::base::GridType ModLagrangeNotAKnotSplineGrid::getType() {
  return sgpp::base::GridType::ModLagrangeNotAKnotSpline;
}

const SBasis& ModLagrangeNotAKnotSplineGrid::getBasis() { return *basis_; }

size_t ModLagrangeNotAKnotSplineGrid::getDegree() { return this->degree; }

Grid* ModLagrangeNotAKnotSplineGrid::unserialize(std::istream& istr) {
  return new ModLagrangeNotAKnotSplineGrid(istr);
}

void ModLagrangeNotAKnotSplineGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModLagrangeNotAKnotSplineGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
