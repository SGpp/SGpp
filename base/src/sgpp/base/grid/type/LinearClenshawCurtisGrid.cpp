// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>

namespace sgpp {
namespace base {

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(std::istream& istr)
    : Grid(istr), generator(storage), boundaryLevel(0) {
  istr >> boundaryLevel;
  generator.setBoundaryLevel(boundaryLevel);
}

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(size_t dim, level_t boundaryLevel)
    : Grid(dim), generator(storage, boundaryLevel), boundaryLevel(boundaryLevel) {}

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(BoundingBox& BB, level_t boundaryLevel)
    : Grid(BB), generator(storage, boundaryLevel), boundaryLevel(boundaryLevel) {}

LinearClenshawCurtisGrid::~LinearClenshawCurtisGrid() {}

sgpp::base::GridType LinearClenshawCurtisGrid::getType() {
  return sgpp::base::GridType::LinearClenshawCurtis;
}

const SBasis& LinearClenshawCurtisGrid::getBasis() {
  static SLinearClenshawCurtisBase basis;
  return basis;
}

std::unique_ptr<Grid> LinearClenshawCurtisGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new LinearClenshawCurtisGrid(istr));
}

void LinearClenshawCurtisGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& LinearClenshawCurtisGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
