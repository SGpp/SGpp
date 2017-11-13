// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/LagrangeSplineBoundaryGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

LagrangeSplineBoundaryGrid::LagrangeSplineBoundaryGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16), boundaryLevel(0) {
  istr >> degree;
  istr >> boundaryLevel;
  basis_.reset(new SLagrangeSplineBase(degree));
  generator.setBoundaryLevel(boundaryLevel);
}

LagrangeSplineBoundaryGrid::LagrangeSplineBoundaryGrid(
    size_t dim, size_t degree, level_t boundaryLevel)
    : Grid(dim),
      generator(storage, boundaryLevel),
      degree(degree),
      basis_(new SLagrangeSplineBase(degree)),
      boundaryLevel(boundaryLevel) {}

LagrangeSplineBoundaryGrid::~LagrangeSplineBoundaryGrid() {}

sgpp::base::GridType LagrangeSplineBoundaryGrid::getType() {
  return sgpp::base::GridType::LagrangeSplineBoundary;
}

const SBasis& LagrangeSplineBoundaryGrid::getBasis() { return *basis_; }

size_t LagrangeSplineBoundaryGrid::getDegree() { return this->degree; }

Grid* LagrangeSplineBoundaryGrid::unserialize(std::istream& istr) {
  return new LagrangeSplineBoundaryGrid(istr);
}

void LagrangeSplineBoundaryGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& LagrangeSplineBoundaryGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
