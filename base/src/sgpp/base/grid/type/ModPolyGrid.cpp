// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

ModPolyGrid::ModPolyGrid(std::istream& istr) : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SPolyModifiedBase(degree));
}

ModPolyGrid::ModPolyGrid(size_t dim, size_t degree)
    : Grid(dim), generator(storage), degree(degree), basis_(new SPolyModifiedBase(degree)) {}

ModPolyGrid::~ModPolyGrid() {}

sgpp::base::GridType ModPolyGrid::getType() { return sgpp::base::GridType::ModPoly; }

SBasis& ModPolyGrid::getBasis() { return *basis_; }

size_t ModPolyGrid::getDegree() const { return this->degree; }

std::unique_ptr<Grid> ModPolyGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new ModPolyGrid(istr));
}

void ModPolyGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModPolyGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
