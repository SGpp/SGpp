// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModLinearClenshawCurtisGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedClenshawCurtisBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {

ModLinearClenshawCurtisGrid::ModLinearClenshawCurtisGrid(std::istream& istr)
    : Grid(istr), generator(storage) {}

ModLinearClenshawCurtisGrid::ModLinearClenshawCurtisGrid(size_t dim)
    : Grid(dim), generator(storage) {
  std::vector<BoundingBox1D> boundingBox1Ds(dim, BoundingBox1D());
  std::vector<Stretching1D> stretching1Ds(dim, Stretching1D("cc"));
  Stretching stretching(boundingBox1Ds, stretching1Ds);
  storage.setStretching(stretching);
}

ModLinearClenshawCurtisGrid::~ModLinearClenshawCurtisGrid() {}

sgpp::base::GridType ModLinearClenshawCurtisGrid::getType() {
  return sgpp::base::GridType::ModLinearClenshawCurtis;
}

SBasis& ModLinearClenshawCurtisGrid::getBasis() { return basis_; }

Grid* ModLinearClenshawCurtisGrid::unserialize(std::istream& istr) {
  return new ModLinearClenshawCurtisGrid(istr);
}

void ModLinearClenshawCurtisGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModLinearClenshawCurtisGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
