// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

NakBsplineGrid::NakBsplineGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SNakBsplineBase(degree));
}

NakBsplineGrid::NakBsplineGrid(size_t dim, size_t degree)
    : Grid(dim), generator(storage), degree(degree), basis_(new SNakBsplineBase(degree)) {}

NakBsplineGrid::~NakBsplineGrid() {}

sgpp::base::GridType NakBsplineGrid::getType() { return sgpp::base::GridType::NakBspline; }

SBasis& NakBsplineGrid::getBasis() { return *basis_; }

size_t NakBsplineGrid::getDegree() { return this->degree; }

Grid* NakBsplineGrid::unserialize(std::istream& istr) { return new NakBsplineGrid(istr); }

void NakBsplineGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& NakBsplineGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
