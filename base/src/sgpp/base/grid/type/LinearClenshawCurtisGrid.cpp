// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(std::istream& istr)
    : Grid(istr), generator(storage) {}

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(size_t dim) : Grid(dim), generator(storage) {
  std::vector<BoundingBox1D> boundingBox1Ds(dim, BoundingBox1D());
  std::vector<Stretching1D> stretching1Ds(dim, Stretching1D("cc"));
  Stretching stretching(boundingBox1Ds, stretching1Ds);
  storage.setStretching(stretching);
}

LinearClenshawCurtisGrid::~LinearClenshawCurtisGrid() {}

sgpp::base::GridType LinearClenshawCurtisGrid::getType() {
  return sgpp::base::GridType::LinearClenshawCurtis;
}

SBasis& LinearClenshawCurtisGrid::getBasis() { return basis_; }

Grid* LinearClenshawCurtisGrid::unserialize(std::istream& istr) {
  return new LinearClenshawCurtisGrid(istr);
}

void LinearClenshawCurtisGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& LinearClenshawCurtisGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
