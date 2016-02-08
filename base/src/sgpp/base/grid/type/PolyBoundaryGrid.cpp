// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>

#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

namespace SGPP {
namespace base {

PolyBoundaryGrid::PolyBoundaryGrid(std::istream& istr) :
  Grid(istr),
  degree(1 << 16),
  basis_(NULL),
  boundaryLevel(0) {
  istr >> degree;
  istr >> boundaryLevel;
}

PolyBoundaryGrid::PolyBoundaryGrid(size_t dim,
                                   size_t degree,
                                   level_t boundaryLevel) :
  Grid(dim),
  degree(degree),
  basis_(NULL),
  boundaryLevel(boundaryLevel) {
}

PolyBoundaryGrid::~PolyBoundaryGrid() {
  if (basis_ != NULL) {
    delete basis_;
  }
}

const SBasis& PolyBoundaryGrid::getBasis() {
  if (basis_ == NULL) {
    basis_ = new SPolyBoundaryBase(degree);
  }

  return *basis_;
}

SGPP::base::GridType PolyBoundaryGrid::getType() {
  return SGPP::base::GridType::PolyBoundary;
}

size_t PolyBoundaryGrid::getDegree() const {
  return this->degree;
}

Grid* PolyBoundaryGrid::unserialize(std::istream& istr) {
  return new PolyBoundaryGrid(istr);
}

void PolyBoundaryGrid::serialize(std::ostream& ostr) {
  this->Grid::serialize(ostr);
  ostr << degree << std::endl;
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* PolyBoundaryGrid::createGridGenerator() {
  return new BoundaryGridGenerator(this->storage, boundaryLevel);
}

}  // namespace base
}  // namespace SGPP
