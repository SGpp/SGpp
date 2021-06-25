// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/NakPBsplineGrid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

NakPBsplineGrid::NakPBsplineGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SNakPBsplineBase(degree));
}

NakPBsplineGrid::NakPBsplineGrid(size_t dim, size_t degree)
    : Grid(dim), generator(storage), degree(degree), basis_(new SNakPBsplineBase(degree)) {}

NakPBsplineGrid::~NakPBsplineGrid() {}

sgpp::base::GridType NakPBsplineGrid::getType() { return sgpp::base::GridType::NakPBspline; }

SBasis& NakPBsplineGrid::getBasis() { return *basis_; }

size_t NakPBsplineGrid::getDegree() { return this->degree; }

Grid* NakPBsplineGrid::unserialize(std::istream& istr) { return new NakPBsplineGrid(istr); }

void NakPBsplineGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& NakPBsplineGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
