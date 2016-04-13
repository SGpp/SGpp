// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

BsplineClenshawCurtisGrid::BsplineClenshawCurtisGrid(std::istream& istr)
    : Grid(istr), generator(storage, boundaryLevel), degree(1 << 16), boundaryLevel(0) {
  istr >> degree;
  istr >> boundaryLevel;
  basis_.reset(new SBsplineClenshawCurtisBase(degree));
}

BsplineClenshawCurtisGrid::BsplineClenshawCurtisGrid(size_t dim, size_t degree,
                                                     level_t boundaryLevel)
    : Grid(dim),
      generator(storage, boundaryLevel),
      degree(degree),
      basis_(new SBsplineClenshawCurtisBase(degree)),
      boundaryLevel(boundaryLevel) {}

BsplineClenshawCurtisGrid::~BsplineClenshawCurtisGrid() {}

sgpp::base::GridType BsplineClenshawCurtisGrid::getType() {
  return sgpp::base::GridType::BsplineClenshawCurtis;
}

const SBasis& BsplineClenshawCurtisGrid::getBasis() { return *basis_; }

size_t BsplineClenshawCurtisGrid::getDegree() { return this->degree; }

std::unique_ptr<Grid> BsplineClenshawCurtisGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new BsplineClenshawCurtisGrid(istr));
}

void BsplineClenshawCurtisGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << degree << std::endl;
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& BsplineClenshawCurtisGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
