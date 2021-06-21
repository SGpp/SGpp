// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

NakBsplineBoundaryGrid::NakBsplineBoundaryGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16), boundaryLevel(0) {
  istr >> degree;
  istr >> boundaryLevel;
  basis_.reset(new SNakBsplineBoundaryBase(degree));
  generator.setBoundaryLevel(boundaryLevel);
}

NakBsplineBoundaryGrid::NakBsplineBoundaryGrid(size_t dim, size_t degree, level_t boundaryLevel)
    : Grid(dim),
      generator(storage, boundaryLevel),
      degree(degree),
      basis_(new SNakBsplineBoundaryBase(degree)),
      boundaryLevel(boundaryLevel) {}

NakBsplineBoundaryGrid::~NakBsplineBoundaryGrid() {}

sgpp::base::GridType NakBsplineBoundaryGrid::getType() {
  return sgpp::base::GridType::NakBsplineBoundary;
}

SBasis& NakBsplineBoundaryGrid::getBasis() { return *basis_; }

size_t NakBsplineBoundaryGrid::getDegree() { return this->degree; }

Grid* NakBsplineBoundaryGrid::unserialize(std::istream& istr) {
  return new NakBsplineBoundaryGrid(istr);
}

void NakBsplineBoundaryGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& NakBsplineBoundaryGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
