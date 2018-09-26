// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "NakBsplineModifiedGrid.hpp"

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

NakBsplineModifiedGrid::NakBsplineModifiedGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SBsplineModifiedBase(degree));
}

NakBsplineModifiedGrid::NakBsplineModifiedGrid(size_t dim, size_t degree)
    : Grid(dim), generator(storage), degree(degree), basis_(new SBsplineModifiedBase(degree)) {}

NakBsplineModifiedGrid::~NakBsplineModifiedGrid() {}

sgpp::base::GridType NakBsplineModifiedGrid::getType() {
  return sgpp::base::GridType::NakBsplineModified;
}

SBasis& NakBsplineModifiedGrid::getBasis() { return *basis_; }

size_t NakBsplineModifiedGrid::getDegree() { return this->degree; }

Grid* NakBsplineModifiedGrid::unserialize(std::istream& istr) {
  return new NakBsplineModifiedGrid(istr);
}

void NakBsplineModifiedGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& NakBsplineModifiedGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
