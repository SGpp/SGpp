// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModNotAKnotBsplineGrid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

ModNotAKnotBsplineGrid::ModNotAKnotBsplineGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SNotAKnotBsplineModifiedBase(degree));
}

ModNotAKnotBsplineGrid::ModNotAKnotBsplineGrid(size_t dim, size_t degree)
    : Grid(dim), generator(storage), degree(degree),
      basis_(new SNotAKnotBsplineModifiedBase(degree)) {}

ModNotAKnotBsplineGrid::~ModNotAKnotBsplineGrid() {}

sgpp::base::GridType ModNotAKnotBsplineGrid::getType() {
  return sgpp::base::GridType::ModNotAKnotBspline;
}

const SBasis& ModNotAKnotBsplineGrid::getBasis() { return *basis_; }

size_t ModNotAKnotBsplineGrid::getDegree() { return this->degree; }

Grid* ModNotAKnotBsplineGrid::unserialize(std::istream& istr) {
  return new ModNotAKnotBsplineGrid(istr);
}

void ModNotAKnotBsplineGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModNotAKnotBsplineGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
