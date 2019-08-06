// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisBoundaryGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>
#include <vector>

namespace sgpp {
namespace base {

LinearClenshawCurtisBoundaryGrid::LinearClenshawCurtisBoundaryGrid(std::istream& istr)
    : Grid(istr), generator(storage), boundaryLevel(0) {
  istr >> boundaryLevel;
  generator.setBoundaryLevel(boundaryLevel);
}

LinearClenshawCurtisBoundaryGrid::LinearClenshawCurtisBoundaryGrid(size_t dim,
                                                                  level_t boundaryLevel)
    : Grid(dim), generator(storage, boundaryLevel), boundaryLevel(boundaryLevel) {
  std::vector<BoundingBox1D> boundingBox1Ds(dim, BoundingBox1D());
  std::vector<Stretching1D> stretching1Ds(dim, Stretching1D("cc"));
  Stretching stretching(boundingBox1Ds, stretching1Ds);
  storage.setStretching(stretching);
}

LinearClenshawCurtisBoundaryGrid::~LinearClenshawCurtisBoundaryGrid() {}

sgpp::base::GridType LinearClenshawCurtisBoundaryGrid::getType() {
  return sgpp::base::GridType::LinearClenshawCurtisBoundary;
}

SBasis& LinearClenshawCurtisBoundaryGrid::getBasis() {
  static SLinearClenshawCurtisBoundaryBase basis;
  return basis;
}

Grid* LinearClenshawCurtisBoundaryGrid::unserialize(std::istream& istr) {
  return new LinearClenshawCurtisBoundaryGrid(istr);
}

void LinearClenshawCurtisBoundaryGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
  ostr << boundaryLevel << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& LinearClenshawCurtisBoundaryGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
