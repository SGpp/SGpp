// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

ModBsplineGrid::ModBsplineGrid(std::istream& istr)
    : Grid(istr), generator(storage), degree(1 << 16) {
  istr >> degree;
  basis_.reset(new SBsplineModifiedBase(degree));
}

ModBsplineGrid::ModBsplineGrid(size_t dim, size_t degree)
    : Grid(dim), generator(storage), degree(degree), basis_(new SBsplineModifiedBase(degree)) {}

ModBsplineGrid::~ModBsplineGrid() {}

sgpp::base::GridType ModBsplineGrid::getType() { return sgpp::base::GridType::ModBspline; }

const SBasis& ModBsplineGrid::getBasis() { return *basis_; }

size_t ModBsplineGrid::getDegree() { return this->degree; }

std::unique_ptr<Grid> ModBsplineGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new ModBsplineGrid(istr));
}

void ModBsplineGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModBsplineGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
