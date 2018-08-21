// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

namespace sgpp {
namespace base {

PolyBoundaryGrid::PolyBoundaryGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16), boundaryLevel(0) {
  istr >> degree;
  istr >> boundaryLevel;
  basis_.reset(new SPolyBoundaryBase(degree));
  generator.setBoundaryLevel(boundaryLevel);
}

PolyBoundaryGrid::PolyBoundaryGrid(size_t dim, size_t degree, level_t boundaryLevel)
    : Grid(dim),
      generator(storage, boundaryLevel),
      degree(degree),
      basis_(new SPolyBoundaryBase(degree)),
      boundaryLevel(boundaryLevel) {}

PolyBoundaryGrid::~PolyBoundaryGrid() {}

SBasis& PolyBoundaryGrid::getBasis() { return *basis_; }

sgpp::base::GridType PolyBoundaryGrid::getType() { return sgpp::base::GridType::PolyBoundary; }

size_t PolyBoundaryGrid::getDegree() const { return this->degree; }

Grid* PolyBoundaryGrid::unserialize(std::istream& istr) {
  return new PolyBoundaryGrid(istr);
}

void PolyBoundaryGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& PolyBoundaryGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
