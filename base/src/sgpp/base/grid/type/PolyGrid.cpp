// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

PolyGrid::PolyGrid(std::istream& istr) :
  Grid(istr),
  generator(storage),
  degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SPolyBase(degree));
}

PolyGrid::PolyGrid(size_t dim, size_t degree) :
  Grid(dim),
  generator(storage),
  degree(degree),
  basis_(new SPolyBase(degree)) {
}

PolyGrid::~PolyGrid() {
}

sgpp::base::GridType PolyGrid::getType() {
  return sgpp::base::GridType::Poly;
}

const SBasis& PolyGrid::getBasis() {
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
GridGenerator& PolyGrid::getGenerator() {
  return generator;
}

}  // namespace base
}  // namespace sgpp
