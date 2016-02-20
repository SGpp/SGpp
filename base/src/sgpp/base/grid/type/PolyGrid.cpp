// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

PolyGrid::PolyGrid(std::istream& istr) :
  Grid(istr),
  degree(1 << 16),
  basis_(NULL) {
  istr >> degree;
}

PolyGrid::PolyGrid(size_t dim, size_t degree) :
  Grid(dim),
  degree(degree),
  basis_(NULL) {
}

PolyGrid::~PolyGrid() {
}

SGPP::base::GridType PolyGrid::getType() {
  return SGPP::base::GridType::Poly;
}

const SBasis& PolyGrid::getBasis() {
  if (basis_ == NULL) {
    basis_ = new SPolyBase(degree);
  }

  return *basis_;
}

size_t PolyGrid::getDegree() const {
  return this->degree;
}

std::unique_ptr<Grid> PolyGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new PolyGrid(istr));
}

void PolyGrid::serialize(std::ostream& ostr) {
  this->Grid::serialize(ostr);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* PolyGrid::createGridGenerator() {
  return new StandardGridGenerator(this->storage);
}

}  // namespace base
}  // namespace SGPP
