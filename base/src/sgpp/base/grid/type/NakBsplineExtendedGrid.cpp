// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

NakBsplineExtendedGrid::NakBsplineExtendedGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SNakBsplineExtendedBase(degree));
}

NakBsplineExtendedGrid::NakBsplineExtendedGrid(size_t dim, size_t degree)
    : Grid(dim), generator(storage), degree(degree), basis_(new SNakBsplineExtendedBase(degree)) {}

NakBsplineExtendedGrid::~NakBsplineExtendedGrid() {}

sgpp::base::GridType NakBsplineExtendedGrid::getType() {
  return sgpp::base::GridType::NakBsplineExtended;
}

SBasis& NakBsplineExtendedGrid::getBasis() { return *basis_; }

size_t NakBsplineExtendedGrid::getDegree() { return this->degree; }

Grid* NakBsplineExtendedGrid::unserialize(std::istream& istr) {
  return new NakBsplineExtendedGrid(istr);
}

void NakBsplineExtendedGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& NakBsplineExtendedGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
